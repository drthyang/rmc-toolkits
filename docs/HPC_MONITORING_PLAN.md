# Remote / HPC Run Monitoring — Design Plan

**Status:** planning (no code yet). **Audience:** rmc-toolkits developers **and** the RMCProfile
team (for the integration conversation).

## 1. Goal

Run RMCProfile on an HPC cluster and monitor the run **remotely** from this app, with the **same
feature set available today** — dashboard R-values / convergence, Atomic Density (KDE / 3D), PCA
Ellipsoid, the AI assistant, and Live Data auto-refresh — without copying the run down by hand.

### Guiding rules (from the product owner)

1. **Security first, performance second.** Nothing in this feature may weaken HPC security or the
   user's credentials.
2. **Listen-only / read-only.** The app **monitors**; it never pushes data, never writes into the
   run, and cannot perturb or control the simulation. The only bytes we send are the minimal
   *read* commands needed to pull output (e.g. `stat`, `tail`, `rsync --read-only`). We treat this
   as a hard boundary: **no write, no submit, no cancel, no parameter change** in v1.
3. **No third parties.** Run data flows only between the HPC and the user's own machine over a
   channel the user already trusts. No cloud relay, no telemetry, no external broker.
4. **Auth = the user's existing SSH trust.** Public-key auth (via the user's `ssh-agent` /
   `~/.ssh/config`), including bastion/jump hosts and MFA/OTP prompts. The app never sees, stores,
   or transmits passwords or private keys.

## 2. Why "listen to a port" becomes "pull over SSH"

The literal idea — RMCProfile opens a port and the app listens — does not survive HPC reality:

- **Compute nodes are firewalled and ephemeral.** They rarely accept inbound connections, get torn
  down when the job ends, and are reached only through the scheduler. A listening port on the node
  is usually unreachable from the outside and often prohibited by center policy.
- **Login nodes don't let users open arbitrary public listeners** either.
- **The one channel every HPC user already has is SSH** (port 22, outbound from the user's machine,
  key-authenticated, MFA-aware, bastion-aware).

So the robust design is **pull-based monitoring over SSH**: the user's local machine reaches *out*
to the cluster over SSH and reads the files RMCProfile already writes. This honors every rule above,
needs **zero inbound exposure**, and — importantly — **works with RMCProfile as it exists today**.

A true low-latency **push** path is kept as an *optional later phase* (§7), also over SSH, never as
a raw public port.

## 3. Architecture (recommended: Option A — SSH pull via the local backend)

```
┌─────────────────┐     SSH (22), outbound, key/MFA auth      ┌─────────────────────┐
│  User's machine │  ───────────────────────────────────────▶ │  HPC login node     │
│                 │      read-only: stat / tail / rsync        │  (run directory on  │
│  rmc-toolkits   │ ◀─────────────────────────────────────── │   shared filesystem)│
│  Flask backend  │      only appended/changed bytes           └────────────────────┘
│  (localhost)    │
│      │          │   local delta cache (…/remote-cache/<run>/)
│      ▼          │
│  existing APIs  │  /api/structure, /api/files, KDE, PCA … unchanged, over the cache
│      │          │
│      ▼          │
│  React frontend │  Dashboard / Atomic Density / PCA Ellipsoid / AI assistant
└─────────────────┘
```

- The **optional Flask backend** (already part of this repo, runs on the user's laptop/workstation,
  binds to localhost only) grows a **remote transport module**. It opens an SSH connection to the
  HPC login node and incrementally syncs the run directory into a **local cache directory**.
- Every existing feature then reads that cache exactly like a local run folder — **no changes to the
  dashboard, KDE, PCA, or assistant code paths.** Heavy compute (KDE / PCA) stays local, as today;
  only raw output files cross the wire.
- **Live Data** maps directly onto adaptive polling of the remote directory.

### Transport choices (all over the same SSH/key infra)

| Transport | Use | Notes |
|---|---|---|
| `rsync --read-only -z` over SSH | whole-directory delta sync | best default; transfers only changed blocks, compresses |
| `sftp` / `ssh cat`,`tail -c +N` | incremental log/CSV tailing | fetch only appended bytes since last offset |
| `sshfs` mount (optional) | reuse the *local folder* code path verbatim | convenient; depends on FUSE availability |

### Why not the alternatives

- **Cloud relay / WebSocket broker:** rejected — data would transit a third party and add an
  inbound attack surface. Violates rules 1 and 3.
- **Reverse tunnel from the compute node (`ssh -R`):** deferred to §7 as an *opt-in* enhancement;
  many centers disallow it and it widens the security surface, so it is never the default.

## 4. Security model (must-haves)

- **Delegate all auth to the system SSH stack.** Shell out to `ssh` / `rsync` / `sftp` (or a library
  that talks to `ssh-agent`), so passphrase-protected keys, `~/.ssh/config`, `ProxyJump` bastions,
  and MFA/OTP prompts all "just work." **The app stores no secrets** — no passwords, no private
  keys, no tokens. (Consistent with this project's rule of never handling credentials in-app.)
- **Read-only by construction.** The remote command set is a fixed allowlist of read operations
  (`stat`, `ls`, `cat`, `tail`, `rsync` in read-only/pull mode). No write, exec-of-arbitrary, or
  scheduler-mutating commands are ever issued.
- **Path allowlist.** Only the user-specified remote run directory (and below) is readable — the
  same data-root restriction the local backend already enforces, applied to the remote side.
- **Strict host-key verification.** Honor `known_hosts`; never auto-accept unknown host keys —
  surface an explicit trust prompt in the UI instead.
- **No inbound exposure.** The backend keeps binding to localhost; SSH is strictly outbound. No new
  listeners, no port-forwarding-in.
- **Least data at rest.** The local cache lives under a user-controlled directory and is clearable;
  document that monitored data is a local copy.

## 5. Performance model (second priority, still important)

- **Connection reuse:** SSH `ControlMaster` + `ControlPersist` so repeated polls share one
  authenticated connection — no re-handshake, and MFA is prompted at most once per session.
- **Incremental transfer:** keep per-file `(size, mtime, head-hash)`; pull **only** appended log
  bytes (`tail -c +offset`) and re-pull large files (e.g. `.rmc6f`, which can be tens of MB for
  ~50k-atom supercells) **only when mtime changes**.
- **Compression:** `rsync -z` / `ssh -C` for the big structure files.
- **Adaptive polling:** poll faster during active convergence, back off when idle; align cadence
  with RMCProfile's write frequency (see §6.2). A user-facing "structure updates" throttle lets
  slow-link users decouple the heavy `.rmc6f`/3D refresh from the light R-value refresh.
- **Compute stays local:** the wire only ever carries raw RMCProfile output; KDE/PCA run on the
  user's machine as they do now.

## 6. Work split

### 6.1 rmc-toolkits (our side) — needed for the MVP

The MVP needs **no RMCProfile changes** — it tails the files RMCProfile already writes.

1. **Backend remote transport** (`web_app/backend/remote/`): SSH/rsync/sftp client, `ssh-agent`
   based, `ControlMaster` reuse, incremental delta sync into a local cache; strict host-key checks;
   read-only command allowlist; remote path allowlist.
2. **Backend API**: `POST /api/remote/connect` (validate host/key, confirm host key, resolve run
   dir), `POST /api/remote/sync` (delta pull, return manifest), `GET /api/remote/status`
   (connected / last-sync / latency). Reuse existing `/api/structure`, `/api/files`, KDE, and PCA
   endpoints unchanged against the cache dir. All localhost-only.
3. **Frontend**: a **"Remote" run source** beside *Select Folder* / *Demo* — a small connection
   form (host, user, remote path, optional jump host) mapped to an SSH profile — plus a status chip
   (connected · last sync · latency) and a persistent **read-only** banner. Wire it into the
   existing run-source abstraction so downstream pages are untouched. Reuse Live Data semantics.
4. **Security hardening + docs**: threat model, "no credentials stored," known-hosts trust UX,
   clear "this is a local copy" messaging.
5. **Tests**: fixture / mock-SSH-server tests for delta sync, offset tailing, host-key handling,
   and read-only enforcement.
6. **Document the static-mode limitation**: the browser-only (GitHub Pages) build **cannot** open
   SSH (no raw sockets in the browser). Remote monitoring therefore **requires the local Flask
   backend**. A browser-only path would need a server-side bridge — out of scope and against rule 3.

### 6.2 RMCProfile team — enhancements (none strictly required for the file-tailing MVP)

These make monitoring cheaper, richer, and more robust. Ordered by value:

1. **Machine-readable status / heartbeat file (highest value).** A small JSON file the run updates
   **atomically** every *N* steps, e.g.:
   ```json
   {
     "schema": 1,
     "run_id": "gts_250k_run3",
     "phase": "running",             // queued | running | finished | failed
     "step": 1250000,
     "moves": {"attempted": 1250000, "accepted": 812004},
     "chi2": {"total": 3.42, "per_dataset": {"neutron": 2.10, "xray": 1.32}},
     "rwp": {"total": 5.8},
     "temperature": 250.0,
     "wall_seconds": 88231,
     "eta_seconds": 42100,
     "updated_unix": 1752570000
   }
   ```
   This lets convergence monitoring read one tiny file instead of scraping several logs, and makes
   "queued / running / finished / failed" unambiguous.
2. **Atomic, consistent writes.** Write partial outputs to a temp file and `rename()` into place (or
   document a flush interval), so a monitor never reads a half-written `.rmc6f` / CSV.
3. **Documented, versioned output schema.** Stable, versioned formats + a documented list of output
   files and their meaning, so remote parsers don't break across RMCProfile releases. (A `schema`
   field as above is the minimal version marker.)
4. **Configurable output directory & write cadence** via the command file, so users can point the
   monitor at one place and trade write frequency against run performance.
5. **(Optional, later) Job-side push agent.** An opt-in mode where RMCProfile — or a companion
   script — streams the compact status feed to a FIFO/socket that the user's local agent tails over
   an **SSH reverse tunnel**, for near-real-time updates without polling. Must be opt-in,
   key-authenticated, and read-only from the app's side. See §7.
6. **(Optional) Scheduler mapping.** A documented job-ID ↔ run-directory convention so the app can
   also show SLURM/PBS status (`squeue` / `sacct`, read-only over the same SSH channel).

## 7. Phasing

- **Phase 8a — MVP (our side only, zero RMCProfile changes).** SSH pull of existing output files →
  the full current feature set, read-only, over the user's SSH trust. Ships independently.
- **Phase 8b — status file.** Consume the RMCProfile JSON heartbeat (§6.2.1) for cheaper, more
  robust convergence and phase reporting.
- **Phase 8c — low-latency / scheduler (optional).** Opt-in reverse-tunnel push for near-real-time
  updates; scheduler status via `squeue`/`sacct`.

## 8. Open questions for the RMCProfile conversation

1. Which HPC auth patterns should we support first — plain SSH key, `ProxyJump` bastion, MFA/OTP,
   Kerberos/GSI (`gsissh`)?
2. Typical output **write cadence** and **file sizes** at production scale; is atomic write feasible?
3. Willingness to emit the **JSON status/heartbeat file**, and can we agree the schema in §6.2.1?
4. Any objection to an **opt-in reverse-tunnel push** later, and what are center policies on it?
5. Output-format **versioning** so parsers stay stable across RMCProfile releases.
6. Is there an existing job-ID ↔ run-directory convention we can lean on for scheduler status?

## 9. One-line summary for the RMCProfile team

> We can monitor RMCProfile HPC runs **today** by pulling the output files read-only over the user's
> own SSH connection — no ports opened, nothing written back, no third parties. The single most
> valuable thing RMCProfile could add is a small, atomically-updated **JSON status/heartbeat file**
> (step, χ²/Rwp, phase, ETA); everything else is an optional performance nicety.

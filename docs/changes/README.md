# Change fragments

Each pull request adds one file to `unreleased/`, named `<short-slug>.md`, holding its changelog
entry exactly as it will appear under `## Unreleased` in [CHANGELOG.md](../CHANGELOG.md): a `###`
heading and its bullets. PRs do not edit CHANGELOG.md itself, so parallel PRs never conflict there.
Before a release, one "fold" PR moves the fragments into the changelog's unreleased section, in merge
order, and deletes them.

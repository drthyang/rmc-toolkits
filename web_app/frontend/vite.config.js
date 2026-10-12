/* global process */
import { defineConfig, searchForWorkspaceRoot } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  base: process.env.VITE_BASE_PATH || '/',
  plugins: [react()],
  // The AI Copilot's skills live at the repository root (.claude/skills),
  // shared with Claude Code; the dev server may read them there.
  server: {
    fs: {
      allow: [searchForWorkspaceRoot(process.cwd()), '../../.claude/skills'],
    },
  },
  test: {
    setupFiles: ['./vitest.setup.js'],
  },
})

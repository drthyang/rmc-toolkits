import js from '@eslint/js'
import globals from 'globals'
import reactHooks from 'eslint-plugin-react-hooks'
import reactRefresh from 'eslint-plugin-react-refresh'
import { defineConfig, globalIgnores } from 'eslint/config'

export default defineConfig([
  globalIgnores(['dist']),
  {
    files: ['**/*.{js,jsx}'],
    extends: [
      js.configs.recommended,
      reactHooks.configs.flat.recommended,
      reactRefresh.configs.vite,
    ],
    languageOptions: {
      ecmaVersion: 2020,
      globals: globals.browser,
      parserOptions: {
        ecmaVersion: 'latest',
        ecmaFeatures: { jsx: true },
        sourceType: 'module',
      },
    },
    rules: {
      'no-unused-vars': ['error', { varsIgnorePattern: '^[A-Z_]', ignoreRestSiblings: true }],
    },
  },
  {
    // The validation harness runs under Node (vite-node).
    files: ['validation/**/*.js'],
    languageOptions: { globals: { ...globals.node } },
  },
  {
    // The run module's in-process adapter (src/openRunFolder.js) is for Node and
    // tests only: the worker modules it loads install self.onmessage whenever
    // self.postMessage is a function, which a window's is. No app module may
    // import it; tests and validation/ may.
    files: ['src/**/*.{js,jsx}'],
    ignores: ['src/**/__tests__/**', 'src/**/*.test.{js,jsx}'],
    rules: {
      'no-restricted-imports': ['error', {
        patterns: [{
          group: ['**/openRunFolder', '**/openRunFolder.js'],
          message: 'openRunFolder is the run module\'s in-process adapter (Node and tests only); the app uses useCopilotRun.',
        }],
      }],
      'no-restricted-syntax': ['error', {
        selector: 'ImportExpression[source.value=/openRunFolder(\\.js)?$/]',
        message: 'openRunFolder is the run module\'s in-process adapter (Node and tests only); the app uses useCopilotRun.',
      }],
    },
  },
])

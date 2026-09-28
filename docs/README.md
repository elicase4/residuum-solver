# Documentation

This directory holds documentation for humans reading top-to-bottom; `CLAUDE.md` at the repo root
stays the terse structural reference (library layering, naming conventions, file-by-file lookup).
When the two overlap, `CLAUDE.md` is the source of truth for *what the convention is* and
`docs/design/` is where *why* and *how to extend it* live.

- **`design/`** — architecture walkthroughs and step-by-step guides for extending the codebase
  (adding a new equation, adding a new solver method, the benchmarking module). Start with
  `CONTRIBUTING.md` at the repo root, which links into this directory.
- **`examples/`** — a catalog of every config under `examples/`: what it demonstrates, what parts
  of the codebase it exercises.
- **`math/`** — PDE derivations (strong form → weak form → discretization) per equation, in
  Markdown or embedded LaTeX. These are intentionally left as stubs for now.

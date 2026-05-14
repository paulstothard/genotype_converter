# Repository Notes

This file is intentionally short. The active documentation is split by audience:

- `README.md`: user-facing overview, installation, common workflows, and
  end-to-end examples.
- `docs/genotype-formats.md`: accepted genotype input layouts, conversion
  behavior, and unconvertible-marker policies.
- `docs/database.md`: SQLite database workflow, source-folder layout, source
  maintenance, mixed-manifest conversion, site-only VCF export, and reference
  downloads.
- `docs/development-plan.md`: maintainer-facing status and remaining work. This
  is tracked in the repository and is accessible from GitHub, but it is not a
  release promise.
- `AGENTS.md`: operating notes for AI coding agents, including validation safety
  and implementation details that should not be repeated in user-facing docs.

Avoid adding new planning material here. Put user instructions in the README or
`docs/`, and put future-work/status notes in `docs/development-plan.md`.

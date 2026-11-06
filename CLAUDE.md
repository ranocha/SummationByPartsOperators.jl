# CLAUDE.md

Cross-tool agent instructions live in `AGENTS.md` (single source of truth) and
are imported below. Add only genuinely Claude-Code-specific notes here, above
the import.

## Claude Code: MCP tools may be deferred

In Claude Code, MCP tools may be *deferred*: a system reminder lists only their
names (e.g., `mcp__julia__julia_eval`), and their schemas must be loaded with
`ToolSearch` (e.g., `select:mcp__julia__julia_eval,mcp__julia__julia_restart`)
before the first call. Deferred tools count as available for the rule to run
Julia code in a persistent session (see [Testing](AGENTS.md#testing)). Before
running any Julia code, check the tool list and the system reminders for
`mcp__julia__*` tools. Use the Bash tool for Julia only if they are not usable
or for steps that AGENTS.md assigns to the shell (the `run_agents` setup
script).

@AGENTS.md

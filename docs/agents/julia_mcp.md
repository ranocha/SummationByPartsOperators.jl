# julia-mcp sessions

Read this before the first call to a julia-mcp tool (`julia_eval`,
`julia_restart`, `julia_list_sessions`; in Claude Code, e.g.,
`mcp__julia__julia_eval`). It complements the section "Testing" in
[`AGENTS.md`](../../AGENTS.md).

- Run `julia utils/setup_run_agents.jl` in the shell first, even if the
  description of an MCP tool tells you to never run Julia via the command line.
- Always pass the **absolute path** of `run_agents/` as `env_path`; do not call
  `Pkg.activate`. The working directory of the session is `env_path`, not the
  repository root, so use absolute paths, e.g.,
  `cd("/path/to/SummationByPartsOperators.jl/test") do ... end` for tests.
- `julia_eval` returns only what is printed to `stdout` and `stderr`, not the
  value of the last expression; use `println` or `display`.
- The default `timeout` is 60 seconds, and a call exceeding its timeout **kills
  the session**. Pass a generous `timeout` (e.g., `1800`) for the first call
  and for the first run of a test file, which is dominated by compilation
  (e.g., about 40 seconds for `test/upwind_operators_test.jl`).
- The developer chooses the Julia version and flags. Check `VERSION`,
  `Threads.nthreads()`, and `Base.JLOptions().check_bounds` if they matter
  (e.g., for `ThreadedMode()` or timings). Other flags can be passed via
  `julia_cmd` (e.g., `julia_cmd = "julia --threads=4"`), which replaces the
  session; pass the same `julia_cmd` in all later calls, since a call without
  it replaces the session again by one with the default flags.
- Revise.jl is loaded automatically, so most changes to `src/` are picked up
  without restarting. **Exception: the `@generated` stencil kernels.** Revise
  updates the helper functions that build their code (e.g.,
  `interior_stencil_expression` or `convolve_interior_code` in
  `src/SBP_operators.jl`), but the kernels keep the previously generated code,
  even for element types not used before. After changing such helpers, and
  after updating `run_agents/`, call `julia_restart` with the same `env_path`.
- Run the formatter in a separate session with the absolute path of
  `run_agents/formatter` as `env_path`, e.g.,
  `using JuliaFormatter; println(format(["/path/to/SummationByPartsOperators.jl/src/file.jl"]))`.

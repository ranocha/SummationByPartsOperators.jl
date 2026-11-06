# AGENTS.md

Canonical instructions for AI coding agents working on
SummationByPartsOperators.jl. Claude Code reads this file via the `@AGENTS.md`
import in `CLAUDE.md`; keep this file the single source of truth and put only
tool-specific notes into tool-specific files. Task-specific details live in
`docs/agents/`; read a linked file there when this file tells you to.

When editing this file, only add instructions and commands that have actually
been tested in this repository. Instructions that do not work are worse than no
instructions.

SummationByPartsOperators.jl is a Julia library of summation-by-parts (SBP)
operators (finite difference, Fourier, Legendre, and general matrix-based
operators) for provably stable semidiscretizations. Correctness of the numerics
matters more than anything else: many bugs are silent (e.g., a wrong
coefficient that only reduces the order of accuracy or breaks the SBP property
at the boundary), so verify mathematical changes.

Before editing, inspect `git status` and preserve all existing work. Do not
discard, reset, or stash unrelated changes. The developer often keeps untracked
files (e.g., PDFs of papers, local benchmarks) in the repository root; leave
them alone and never commit them.

## Repository layout and architecture

- `src/SummationByPartsOperators.jl` defines the abstract type hierarchy
  (`AbstractDerivativeOperator{T}` and its subtypes; `T` is the real type,
  e.g., `Float32`, `Float64`, `BigFloat`, or `Rational{BigInt}`), includes all
  source files, and lists all exports. The common interface (`mul!`,
  `mass_matrix`, `grid`, `integrate`, ...) is documented in
  `src/general_operators.jl`.
- Finite difference operators (`DerivativeOperator` and its periodic, upwind,
  dissipation, and variable-coefficient variants) get their coefficients from
  a `source` (subtype of `SourceOfCoefficients`, named after the publication,
  e.g., `MattssonNordström2004()`) via functions such as
  `first_derivative_coefficients(source, order, T, mode)`. Each source lives
  in its own file in `src/SBP_coefficients/` with a docstring and a
  `Base.show` method citing the publication.
- Constructors take an execution mode: `SafeMode()` (plain Julia loops),
  `FastMode()` (default, LoopVectorization.jl `@turbo`), or `ThreadedMode()`
  (`@tturbo`). The stencil kernels are `@generated` functions whose code is
  built by helper functions (e.g., in `src/SBP_operators.jl`); Revise.jl does
  not pick up changes to these helpers, see
  [`docs/agents/julia_mcp.md`](docs/agents/julia_mcp.md).
- `ext/`: package extensions for weak dependencies; `docs/`: Documenter.jl
  documentation (`docs/src/`) and agent instructions (`docs/agents/`).
- `dev/`: scripts used to derive or check coefficients from the literature;
  not part of the package, not tested, and excluded from formatting.
- `run/` (if present, ignored by git): the human developer's own Julia project.
  Do not modify it; agents use `run_agents/` instead (see below).

## Testing

The test files in `test/` are plain scripts with `@testset`s. `test/runtests.jl`
includes them in the partitions `part1`, `part2`, and `rational`, selected by
the environment variable `SBP_TEST` (default `all`); CI runs all of them for
every pull request. Do not run the full test suite (`Pkg.test()`,
`julia test/runtests.jl`, or a whole partition) unless the user asks for it or
a change affects all operators (e.g., the common interface or the stencil
kernels). On the developer's machine, the partitions take about 4 (`part1`),
6 (`part2`), and 2.5 (`rational`) minutes, almost all of it compilation. Run
only the test files related to your change instead.

**Run Julia code (tests, scratch scripts) in a persistent session via
[julia-mcp](https://github.com/aplavin/julia-mcp) if its tools (`julia_eval`,
`julia_restart`, `julia_list_sessions`) are available to you** (possibly only
listed by name until they are loaded, see `CLAUDE.md` for Claude Code). Read
[`docs/agents/julia_mcp.md`](docs/agents/julia_mcp.md) before the first call.
Use the shell commands in this file only if no such session is available.

### The `run_agents` Julia projects

Use the Julia project `run_agents/` (ignored by git) for tests and scratch
scripts. It develops this checkout of SummationByPartsOperators.jl and contains
all test dependencies. The formatter lives in `run_agents/formatter/`. Create
or update both projects **in the shell before running tests or the
formatter**, e.g., from the root of the repository,

```bash
julia utils/setup_run_agents.jl
```

It takes only a few seconds if nothing changed. It re-creates `run_agents/`
(discarding packages you added there) if `test/Project.toml` changed or the
manifest was resolved with another minor version of Julia. If loading the
package fails with `ArgumentError: Package SummationByPartsOperators does not
have X in its dependencies` (after a dependency was added to `Project.toml`),
run `julia --project=run_agents -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'`.
Packages for your own experiments (e.g., BenchmarkTools) may be added to
`run_agents/`, never to `test/Project.toml` or `Project.toml` unless this is
part of the task.

### Running selected tests

Run the Julia code of the following command in the julia-mcp session if
available, otherwise in the shell from the root of the repository:

```bash
julia --project=run_agents -e '
using Test, SummationByPartsOperators
cd("test") do
    @testset "Upwind" include("upwind_operators_test.jl")
    @testset "Legendre" include("legendre_operators_test.jl")
end'
```

- Find the relevant test files by searching `test/` for what you changed,
  e.g., `grep -rn "Mattsson2017" test/`. Many operators are also tested in
  `test/rational_test.jl` (exact arithmetic), `test/special_matrix_types.jl`
  (composite element types), and `test/ad_test.jl` (ForwardDiff.jl).
- Wrap each file in its own `@testset` as above so that a failure does not
  abort the remaining files. Some test files define helper functions at the
  top level; redefining them in the same session is fine.
- `include` test files inside `cd(".../test") do ... end` so that their
  relative `include`s work.
- Do not use `julia --project=.` or `julia --project=test` (each lacks
  packages), and do not stack environments via `JULIA_LOAD_PATH`.

### Adding and changing tests

- Cover new or changed behavior with focused tests in the matching
  `test/*_test.jl` file. For a bug fix, verify that the regression test
  detects the original behavior when feasible (e.g., by temporarily reverting
  only your own changes).
- Follow the existing tests: check the order of accuracy (exactness for
  polynomials at the boundary and in the interior), the SBP property, and the
  agreement of different ways to apply the operator (e.g., `mul!` with and
  without `α`/`β`), usually for `T in (Float32, Float64)`.
- Also test new coefficients that are known as exact rational numbers in
  `test/rational_test.jl` (`Rational{BigInt}`, checked with `==`).
- Never loosen tolerances just to make a failing test pass. A local, commented
  tolerance is acceptable only for small differences whose cause is understood
  and outside the package (e.g., a different LLVM version); tell the user.
- Prefer `@test_throws "message" f(...)` when users rely on the message, in
  addition to or instead of the exception type.

## Code style

- Format all changed Julia files with the JuliaFormatter.jl version pinned in
  `run_agents/formatter/` (v1.0.60, as in CI; other versions format
  differently), e.g., from the root of the repository,
  ```bash
  julia --project=run_agents/formatter -e 'using JuliaFormatter; println(format(["src/file.jl", "test/file_test.jl"]))'
  ```
  `format` returns `true` if the files were already formatted; otherwise, it
  formats them in place. Pass only the files you changed; never run
  `format(".")`, which also formats the developer's untracked files.
- Keep code type-stable and generic in the real type `T` (including
  `BigFloat`, `Rational{BigInt}`, and AD types). Write coefficients known
  exactly as rationals converted to `T` (e.g., `T(-1 // 2)`) so that the SBP
  property holds exactly for `Rational{BigInt}`.
- Respect the axes of input arrays (`eachindex`, `firstindex`, `axes`, ...)
  instead of assuming `1:length(u)`, and check sizes with `@argcheck` or
  informative `DimensionMismatch`/`ArgumentError`s.
- `mul!` must not allocate; check with `@allocated mul!(du, D, u)` after a
  first call for compilation.
- Document and export new public functions and types (the API reference
  includes all docstrings via `@autodocs`). Name new sources of coefficients
  after the publication (authors and year) and cite it in the docstring and
  the `Base.show` method.
- Match the style of the surrounding code, and look for existing helpers
  before adding new ones.

## Pull requests

- Develop in a branch; keep changes focused on a single goal.
- Document user-visible changes, including bug fixes that change results, in
  `NEWS.md` (section "Added", "Changed", or "Fixed" of the current lifecycle).
- Do not change the version number in `Project.toml`; the maintainer does this
  in separate commits.
- Disclose major LLM/AI contributions in the PR description (see
  `CONTRIBUTING.md`).
- Do not commit `Manifest.toml` files, `docs/build/`, PDFs, `run*/`
  directories, or temporary notes.
- Do not commit or push unless asked to do so by the user.

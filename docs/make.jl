using Oscar
local_build = false
bla = normpath(joinpath(dirname(pathof(Oscar)), "..", "docs", "make_work.jl"))
include(bla)

#deploydocs(
#   julia = "1.3",
#   repo   = "github.com/oscar-system/Oscar.jl.git",
#   target = "build",
#   deps = nothing,
#   make   = nothing,
#   osname = "linux"
#)

deploydocs(
   repo   = "github.com/oscar-system/Oscar.jl.git",
#  deps = Deps.pip("pymdown-extensions", "pygments", "mkdocs", "python-markdown-math", "mkdocs-material", "mkdocs-cinder"),
   deps = nothing,
   target = "build",
   push_preview = true,
#  make = () -> run(`mkdocs build`),
   make = nothing
)


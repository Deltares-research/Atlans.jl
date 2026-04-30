using Atlans
using Documenter
using DocumenterMarkdown


makedocs(;
	sitename = "Atlans.jl",
	repo = "https://github.com/Deltares-research/Atlans.jl",
	format = Markdown(),
	authors = "Deltares and contributors",
	modules = [Atlans],
)

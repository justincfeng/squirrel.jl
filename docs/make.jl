using Documenter, squirrel

makedocs(
    sitename    = "squirrel.jl",
    format      = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors     = "Justin C. Feng and Filip Hejda and Sante Carloni",
    pages = [
                "Home"               => "index.md",
                "Metrics"            => "metrics.md",
                "Geodesic solver"   => "geodesics.md",
                "Squirrel algorithm" => "squirrel.md",
                "Evaluation"         => "evaluation.md",
            ]
        )

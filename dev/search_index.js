var documenterSearchIndex = {"docs":
[{"location":"91-developer/#dev_docs","page":"Developer documentation","title":"Developer documentation","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"note: Contributing guidelines\nIf you haven't, please read the Contributing guidelines first.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If you want to make contributions to this package that involves code, then this guide is for you.","category":"page"},{"location":"91-developer/#First-time-clone","page":"Developer documentation","title":"First time clone","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: If you have writing rights\nIf you have writing rights, you don't have to fork. Instead, simply clone and skip ahead. Whenever upstream is mentioned, use origin instead.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If this is the first time you work with this repository, follow the instructions below to clone the repository.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fork this repo\nClone your repo (this will create a git remote called origin)\nAdd this repo as a remote:\ngit remote add upstream https://github.com/mmikhasenko/ThreeBodyDecays.jl","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"This will ensure that you have two remotes in your git: origin and upstream. You will create branches and push to origin, and you will fetch and update your local main branch from upstream.","category":"page"},{"location":"91-developer/#Linting-and-formatting","page":"Developer documentation","title":"Linting and formatting","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Install a plugin on your editor to use EditorConfig. This will ensure that your editor is configured with important formatting settings.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We use https://pre-commit.com to run the linters and formatters. In particular, the Julia code is formatted using JuliaFormatter.jl, so please install it globally first:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # Press ]\npkg> activate\npkg> add JuliaFormatter","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To install pre-commit, we recommend using pipx as follows:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"# Install pipx following the link\npipx install pre-commit","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"With pre-commit installed, activate it as a pre-commit hook:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit install","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To run the linting and formatting manually, enter the command below:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit run -a","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Now, you can only commit if all the pre-commit tests pass.","category":"page"},{"location":"91-developer/#Testing","page":"Developer documentation","title":"Testing","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"As with most Julia packages, you can just open Julia in the repository folder, activate the environment, and run test:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # press ]\npkg> activate .\npkg> test","category":"page"},{"location":"91-developer/#Working-on-a-new-issue","page":"Developer documentation","title":"Working on a new issue","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We try to keep a linear history in this repo, so it is important to keep your branches up-to-date.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fetch from the remote and fast-forward your local main\ngit fetch upstream\ngit switch main\ngit merge --ff-only upstream/main\nBranch from main to address the issue (see below for naming)\ngit switch -c 42-add-answer-universe\nPush the new local branch to your personal remote repository\ngit push -u origin 42-add-answer-universe\nCreate a pull request to merge your remote branch into the org main.","category":"page"},{"location":"91-developer/#Branch-naming","page":"Developer documentation","title":"Branch naming","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If there is an associated issue, add the issue number.\nIf there is no associated issue, and the changes are small, add a prefix such as \"typo\", \"hotfix\", \"small-refactor\", according to the type of update.\nIf the changes are not small and there is no associated issue, then create the issue first, so we can properly discuss the changes.\nUse dash separated imperative wording related to the issue (e.g., 14-add-tests, 15-fix-model, 16-remove-obsolete-files).","category":"page"},{"location":"91-developer/#Commit-message","page":"Developer documentation","title":"Commit message","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Use imperative or present tense, for instance: Add feature or Fix bug.\nHave informative titles.\nWhen necessary, add a body with details.\nIf there are breaking changes, add the information to the commit message.","category":"page"},{"location":"91-developer/#Before-creating-a-pull-request","page":"Developer documentation","title":"Before creating a pull request","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: Atomic git commits\nTry to create \"atomic git commits\" (recommended reading: The Utopic Git History).","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Make sure the tests pass.\nMake sure the pre-commit tests pass.\nFetch any main updates from upstream and rebase your branch, if necessary:\ngit fetch upstream\ngit rebase upstream/main BRANCH_NAME\nThen you can open a pull request and work with the reviewer to address any issues.","category":"page"},{"location":"91-developer/#Building-and-viewing-the-documentation-locally","page":"Developer documentation","title":"Building and viewing the documentation locally","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Following the latest suggestions, we recommend using LiveServer to build the documentation. Here is how you do it:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Run julia --project=docs to open Julia in the environment of the docs.\nIf this is the first time building the docs\nPress ] to enter pkg mode\nRun pkg> dev . to use the development version of your package\nPress backspace to leave pkg mode\nRun julia> using LiveServer\nRun julia> servedocs()","category":"page"},{"location":"91-developer/#Making-a-new-release","page":"Developer documentation","title":"Making a new release","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To create a new release, you can follow these simple steps:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Create a branch release-x.y.z\nUpdate version in Project.toml\nUpdate the CHANGELOG.md:\nRename the section \"Unreleased\" to \"[x.y.z] - yyyy-mm-dd\" (i.e., version under brackets, dash, and date in ISO format)\nAdd a new section on top of it named \"Unreleased\"\nAdd a new link in the bottom for version \"x.y.z\"\nChange the \"[unreleased]\" link to use the latest version - end of line, vx.y.z ... HEAD.\nCreate a commit \"Release vx.y.z\", push, create a PR, wait for it to pass, merge the PR.\nGo back to main screen and click on the latest commit (link: https://github.com/mmikhasenko/ThreeBodyDecays.jl/commit/main)\nAt the bottom, write @JuliaRegistrator register","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"After that, you only need to wait and verify:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Wait for the bot to comment (should take < 1m) with a link to a RP to the registry\nFollow the link and wait for a comment on the auto-merge\nThe comment should said all is well and auto-merge should occur shortly\nAfter the merge happens, TagBot will trigger and create a new GitHub tag. Check on https://github.com/mmikhasenko/ThreeBodyDecays.jl/releases\nAfter the release is create, a \"docs\" GitHub action will start for the tag.\nAfter it passes, a deploy action will run.\nAfter that runs, the stable docs should be updated. Check them and look for the version number.","category":"page"},{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [ThreeBodyDecays]","category":"page"},{"location":"95-reference/#ThreeBodyDecays.ThreeBodyDecay-Tuple{Any}","page":"Reference","title":"ThreeBodyDecays.ThreeBodyDecay","text":"ThreeBodyDecay(descriptor)\n\nConstructs a ThreeBodyDecay object using one argument, a descriptor. The descriptor is a list of pairs, names .=> zip(couplings, chains).\n\nExamples\n\nThreeBodyDecay(\"K892\" .=> zip([1.0, -1.0, 0.2im], [chain1, chain2, chain3]))\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.ThreeBodyDecay-Union{Tuple{S}, Tuple{L}, Tuple{T}, Tuple{Vector{T}, Vector{L}, Vector{S}}} where {T<:AbstractDecayChain, L<:Number, S<:AbstractString}","page":"Reference","title":"ThreeBodyDecays.ThreeBodyDecay","text":"ThreeBodyDecay(; chains, couplings, names)\n\nConstructs a ThreeBodyDecay object with the given parameters.\n\nArguments\n\nchains: An array of chains involved in the decay. The length of this array should match the lengths of couplings and names.\ncouplings: An array of coupling constants for each chain in the decay. The length of this array should match the lengths of chains and names.\nnames: An array of names for each chain, or names of resonances in the decay. The length of this array should match the lengths of chains and couplings.\n\nReturns\n\nA ThreeBodyDecay object with the specified chains, couplings, and names.\n\nExamples\n\nThreeBodyDecay(\nchains=[chain1, chain2, chain3],\ncouplings=[1.0, -1.0, 0.2im],\nnames=[\"L1405\", \"L1405\", \"K892\"])\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#Base.vcat-Tuple{Vararg{ThreeBodyDecay}}","page":"Reference","title":"Base.vcat","text":"Base.vcat(models::ThreeBodyDecay...)\n\nConcatenates multiple ThreeBodyDecay objects into a single ThreeBodyDecay. Argument is variable number of ThreeBodyDecay objects.\n\nAn example\n\nextended_model = vcat(model[2], model[2:3], model)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.DecayChainLS-Tuple{}","page":"Reference","title":"ThreeBodyDecays.DecayChainLS","text":"DecayChainLS(; k, # chain is specified by the spectator index k Xlineshape, # lambda function for lineshape jp, # the spin-parity of the resonance, e.g. jp\"1/2-\" Ps, # need parities, e.g. Ps=ThreeBodyParities('+','+','+'; P0='+') tbs) # give three-body-system structure\n\nReturns the decay chain with the smallest LS, ls\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.DecayChainsLS-Tuple{}","page":"Reference","title":"ThreeBodyDecays.DecayChainsLS","text":"DecayChainsLS(; k, # chain is specified by the spectator index k Xlineshape, # lambda function for lineshape jp, # the spin-parity of the resonance, e.g. jp\"1/2-\" Ps, # need parities, e.g. Ps=ThreeBodyParities('+','+','+'; P0='+') tbs) # give three-body-system structure\n\nReturns an array of the decay chains with all possible couplings\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.Invariants-Union{Tuple{@NamedTuple{m1::T, m2::T, m3::T, m0::T}}, Tuple{T}} where T","page":"Reference","title":"ThreeBodyDecays.Invariants","text":"Invariants(ms::MassTuple{T}; σ1, σ2) Invariants(ms::MassTuple{T}; σ1, σ3) Invariants(ms::MassTuple{T}; σ2, σ3)\n\nConstruct a tuple of (σ1, σ2, σ3) from just two invariants and the mass tuple.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.cosθij-Tuple{Any, Any}","page":"Reference","title":"ThreeBodyDecays.cosθij","text":"cosθij(k,σs,msq)\n\nIsobar decay angle for the chain-k, i.e. an angle of between vectors pi and -pk in the (ij) rest frame.\n\nExplicit forms: cosθ23, cosθ31, cosθ12.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.phase-Tuple{Any}","page":"Reference","title":"ThreeBodyDecays.phase","text":"Phase for wigner d-functions for clockwise rotations\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.phase_space_integrand-Tuple{Any, Any}","page":"Reference","title":"ThreeBodyDecays.phase_space_integrand","text":"phase_space_integrand(function_σs, ms; k::Int)\n\nCalculate the phase space integrand for a given function function_σs, and mass tuple ms. The key argument k specifies the mapping: σk->[0,1], zk->[0,1]. It returns an integrand function of x, x ∈ [0,1]x[0,1] domain to pass to a numerical integrator.\n\nArguments\n\nfunction_σs: A function that takes a MandelstamTuple and returns a scalar.\nms: A scalar representing the mass.\nk: An integer representing the mapping index.\n\nUsage\n\nintegrand = phase_space_integrand(function_σs, ms; k)\n\nSee also\n\nx2σs\nprojection_integrand\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.polardalitz2invariants-Tuple{Any, Tuple}","page":"Reference","title":"ThreeBodyDecays.polardalitz2invariants","text":"polardalitz2invariants(θ, expansion_point)\n\nFor given polar angle θ, it returns an (σ1,σ2,σ3) Tuple of polynomials of radius r(θ) around the expansion point. The polynomial works as a function of the r coordinate.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.projection_integrand-Tuple{Any, Any, Any}","page":"Reference","title":"ThreeBodyDecays.projection_integrand","text":"projection_integrand(function_σs, ms, σk; k)\n\nCalculate the projection integrand for a given function function_σs, mass tuple ms, and Mandelstam variable σk, with k specified by a keyword argument. It returns an integrand function of x, x ∈ [0,1] to pass to a numerical integrator.\n\nArguments\n\nfunction_σs: A function that takes a MandelstamTuple and returns a scalar.\nms: A scalar representing the mass.\nσk: A scalar representing the Mandelstam variable.\nk: A scalar representing the momentum transfer (optional).\n\nUsage\n\nplot(4.2, 4.6) do e1\n\tI = Base.Fix1(unpolarized_intensity, model)\n\tintegrand = projection_integrand(I, masses(model), e1^2; k = 3)\n\te1 * quadgk(integrand, 0, 1)[1]\nend\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.unpolarized_intensity-Tuple{Any, Any}","page":"Reference","title":"ThreeBodyDecays.unpolarized_intensity","text":"unpolarized_intensity(model::ThreeBodyDecay, σs; kw...)\n\nComputes squared amplitude summed over spin projections.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.wr","page":"Reference","title":"ThreeBodyDecays.wr","text":"wr(system_a, reference_b, particle_c=0)\n\nCreate a WignerRotation object of the right type based on provided indices. The daughter particles are numbered 1,2,3, the mother particle is 0.\n\nsystem_a tells which isobar is considered, and\nreference_b tell which system is used as a reference.\n\nFor system_a and reference_b the spectator notations are used, i.e. 1 for the system (2,3), 2 for the system (3,1), and 3 for the system (1,2).\n\n\n\n\n\n","category":"function"},{"location":"95-reference/#ThreeBodyDecays.σiofk-Tuple{Any, Any, Any}","page":"Reference","title":"ThreeBodyDecays.σiofk","text":"σiofk(k,z,σj,msq)\n\nComputes invariant σi = (p0 - pi)² from the scattering angle z=cosθij in the rest from of (i,j), given the mass of the system m(i,j)² = σk\n\nExplicit forms: σ3of2, σ1of3, σ2of1.\n\nSee also σjofk(z,σk,msq; k)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#ThreeBodyDecays.σjofk-Tuple{Any, Any, Any}","page":"Reference","title":"ThreeBodyDecays.σjofk","text":"σjofk(z,σi,msq; k::Int)\n\nComputes invariant σj = (p0-pj)² from the scattering angle z=cosθij in the rest from of (i,j), given the mass of the system m(i,j)² = σk\n\nExplicit forms: σ3of1, σ1of2, σ2of3.\n\nSee also σiofk(z,σj,msq; k)\n\n\n\n\n\n","category":"method"},{"location":"90-contributing/#contributing","page":"Contributing guidelines","title":"Contributing guidelines","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"First of all, thanks for the interest!","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"We welcome all kinds of contribution, including, but not limited to code, documentation, examples, configuration, issue creating, etc.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"Be polite and respectful, and follow the code of conduct.","category":"page"},{"location":"90-contributing/#Bug-reports-and-discussions","page":"Contributing guidelines","title":"Bug reports and discussions","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you think you found a bug, feel free to open an issue. Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.","category":"page"},{"location":"90-contributing/#Working-on-an-issue","page":"Contributing guidelines","title":"Working on an issue","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you found an issue that interests you, comment on that issue what your plans are. If the solution to the issue is clear, you can immediately create a pull request (see below). Otherwise, say what your proposed solution is and wait for a discussion around it.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"tip: Tip\nFeel free to ping us after a few days if there are no responses.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If your solution involves code (or something that requires running the package locally), check the developer documentation. Otherwise, you can use the GitHub interface directly to create your pull request.","category":"page"},{"location":"","page":"ThreeBodyDecays","title":"ThreeBodyDecays","text":"CurrentModule = ThreeBodyDecays","category":"page"},{"location":"#ThreeBodyDecays","page":"ThreeBodyDecays","title":"ThreeBodyDecays","text":"","category":"section"},{"location":"","page":"ThreeBodyDecays","title":"ThreeBodyDecays","text":"Documentation for ThreeBodyDecays.","category":"page"},{"location":"#Contributors","page":"ThreeBodyDecays","title":"Contributors","text":"","category":"section"},{"location":"","page":"ThreeBodyDecays","title":"ThreeBodyDecays","text":"<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->\n<!-- prettier-ignore-start -->\n<!-- markdownlint-disable -->\n\n<!-- markdownlint-restore -->\n<!-- prettier-ignore-end -->\n\n<!-- ALL-CONTRIBUTORS-LIST:END -->","category":"page"}]
}
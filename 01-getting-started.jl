### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2d99f63e-5489-11ed-33b9-4d6021e041dd
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("ImageMorphology")
		Pkg.add("GLM")
		Pkg.add("Statistics")
		Pkg.add(url="https://github.com/Dale-Black/DICOM.jl")
		Pkg.add(url="https://github.com/Dale-Black/DICOMUtils.jl")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI, CSV, DataFrames, CairoMakie, ImageMorphology, GLM, Statistics, DICOM, DICOMUtils, CalciumScoring
end

# ╔═╡ c0bf04c7-7655-469a-8aaf-68f278bd75f6
md"""
# Overview
In this tutorial, we will briefly showcase the different calcium scoring techniques available in CalciumScoring.jl. Specifically **(1) Agatston Scoring**, **(2) Integrated Calcium Scoring**, and **(3) Spatially Weighted Calcium Scoring**
"""

# ╔═╡ cc4ffe96-d49e-4ed2-955a-14b7d0617757
md"""
## Import Packages
First, let's import the most up-to-date version of CalciumScoring.jl, which can be found on the main/master branch of the [GitHub repository](https://github.com/Dale-Black/CalciumScoring.jl).

Because we are using the unregistered version (most recent), we will need to `Pkg.add` this explicitly without using Pluto's built-in package manager. Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.

To help with the formatting of this documentation, we will also add [PlutoUI.jl](https://github.com/JuliaPluto/PlutoUI.jl). We will also use a local fork of [DICOM.jl]() along with a simple package for working with DICOM images, DICOMUtils.jl. Lastly, to visualize the results, we're going to add [Makie.jl](https://github.com/JuliaPlots/Makie.jl).
"""

# ╔═╡ ad36247d-9f0c-4083-bdee-21409b647c97
TableOfContents()

# ╔═╡ 830b0add-f62e-4340-8ee9-112f13e99326
md"""
## Load Simulated Phantom Images
To understand calcium scoring, we will load simulated phantom images with calcium inserts of various densities and sizes.
"""

# ╔═╡ 7c0763ec-cb65-40ae-9acf-339995c54284
root_dir = joinpath(pwd(), "images")

# ╔═╡ 9525fdf5-cff1-41d2-b7b3-52f34df6af5a
dicom_dir = joinpath(root_dir, "120")

# ╔═╡ e4c99f57-f128-4759-bf58-02f5845d7d11
mask_dir = joinpath(root_dir, "calcium-insert-masks")

# ╔═╡ 45a05f47-f229-4d97-b5ae-929a27ac589e
begin
	dcm = dcmdir_parse(dicom_dir)
	dcm_array = load_dcm_array(dcm)
	dcm_array_inserts = dcm_array[:, :, 4:end]
end;

# ╔═╡ e1627af9-6dd7-4713-a9d1-e94cc7c4a890
begin
	mask_L_HD = Array(CSV.read(string(mask_dir, "/mask_L_HD.csv"), DataFrame; header=false))
end;

# ╔═╡ decb80cf-e696-47fd-9d50-04e0579b7fe6
md"""
## Visualize
Below, we can visualize the simulated coronary artery calcification phantom, along with the generated masks (overlayed in red)

Let's first dilate the mask to make sure we include all potential calcium
"""

# ╔═╡ f4b2bd5f-dc8d-45ff-9629-a712c489ee60
begin
    mask_L_HD_3D = Array{Bool}(undef, size(dcm_array_inserts))
    for z in 1:size(dcm_array_inserts, 3)
		mask_L_HD_3D[:, :, z] = mask_L_HD
    end
end;

# ╔═╡ 937bbe02-09c7-4c1e-8770-d5c681784958
dilated_mask_L_HD_3D = dilate(dilate(mask_L_HD_3D));

# ╔═╡ 0afc60d7-9608-4189-b85a-7dc97a498e98
md"""
#### Helper functions
"""

# ╔═╡ 6e6704d8-db30-4a17-ad08-c1928f5c1750
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 83773d29-715f-45d4-a224-57d227c63f4a
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=2, show_value=true)
end

# ╔═╡ 8cb8f78b-515c-46ee-9edc-ea22b538a5b5
function overlay_mask_plot(array, mask, var, title::AbstractString)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    indices_lbl = findall(x -> x == zs[var], label_array[:, 3])

    fig = Figure()
    ax = Makie.Axis(fig[1, 1])
    ax.title = title
    heatmap!(array[:, :, zs[var]]; colormap=:grays)
    scatter!(
        label_array[:, 1][indices_lbl],
        label_array[:, 2][indices_lbl];
        markersize=1,
        color=:red,
    )
    return fig
end

# ╔═╡ bb768658-c623-4941-b248-30fb6fe76779
@bind a overlay_mask_bind(dcm_array_inserts)

# ╔═╡ 6ad2569c-2e98-4abe-a02e-80943e8b7b2f
overlay_mask_plot(permutedims(dcm_array_inserts, (2, 1, 3)), permutedims(dilated_mask_L_HD_3D, (2, 1, 3)), a, "Overlayed Mask")

# ╔═╡ b8738edd-3be2-44fc-9812-a96332c51f5f
md"""
# 1. Agatston Scoring
Calcium scoring is a technique for measuring calcium in the coronary arteries within a CT scan. Calcium scoring is traditionally associated with the Agatston scoring method [(1)](https://pubmed.ncbi.nlm.nih.gov/2407762/)

This package allows users to compute the traditional Agatston score seen below.
"""

# ╔═╡ f07faa8a-6653-46ac-93b1-8e8d6ef53a94
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ 89b45727-74ab-4c3b-bf24-53642def397f
spacing = get_pixel_size(dcm[1].meta)

# ╔═╡ fca9ac32-e940-4a6a-8f84-51547c5329ce
overlayed_mask_l_hd = create_mask(dcm_array_inserts, dilated_mask_L_HD_3D);

# ╔═╡ e27ed808-20cf-4846-b2d3-53ef3dd9ff71
begin
	alg = Agatston()
	score(overlayed_mask_l_hd, spacing, alg)
end

# ╔═╡ d453b50f-06e5-4cab-a6b7-a0f1629770a9
md"""
We can also compute the calcium mass score, via the Agatston technique by computing the mass calibration factor. This factor translates an Agatston score to a mass score via [CITE]()

We could use PhantomSegmentation.jl for this, but an already known estimate of 	0.00075 should work
"""

# ╔═╡ 6bad9fee-b5b1-4c58-a77c-95cff35982ef
md"""
## Results
"""

# ╔═╡ 480dcaca-861a-4d51-8bcc-dc3e24980e9c
begin
	mass_calibration = 0.00075
	agat_score, mass_agat_score = score(overlayed_mask_l_hd, spacing, mass_calibration, alg)
end

# ╔═╡ 1fafecbf-b767-49dc-938f-ad2b6c4241de
md"""
Let's compare that with the ground truth mass of the large, high-density calcium insert
"""

# ╔═╡ 555eac6e-ce5c-4002-8a31-7e39ed0174c3
begin
	density_ground_truth = 800 # mg/cc
	volume_ground_truth = 176.625 # mm^3
	mass_ground_truth = volume_ground_truth * density_ground_truth * 1e-3 # mg
end

# ╔═╡ 00fcc622-cf8b-43ee-a5bd-acafabc903f7
md"""
We see that the ground truth mass (`mass_ground_truth` = $(mass_ground_truth)) is very similar to the calculated mass (`mass_agat_score` = $(mass_agat_score))
"""

# ╔═╡ a4510ff6-92df-4aac-a774-7f046cfa9ea4
md"""
# 2. Integrated Calcium Scoring
For integrated calcium scoring, we need the dilated mask from above, along with a background mask, and a calibration line for ``S_{Obj}``
"""

# ╔═╡ dff1746c-6e43-4cdc-a41d-b3c26ee9e680
md"""
## Calibration Line
For a true calibration line, one would need to include known calcium density inserts in the original image for patient-specific calibration. For this example, we will use previously validated calibration points, specific to these simulated images. This will allow us to then compute a calibration line, using GLM.jl
"""

# ╔═╡ a2d5117a-04e2-4476-9ef6-6d38070a2264
begin
	calibration_densities = [0, 200, 400, 800]
	calibration_intensities = [0, 297.429, 545.245, 1075.82]
	df_calibration = DataFrame(:density => calibration_densities, :intensity => calibration_intensities)
end

# ╔═╡ 682c3d99-11c6-4927-bcf5-0c1833a5588e
linearRegressor = lm(@formula(intensity ~ density), df_calibration);

# ╔═╡ b26f3244-ea93-41aa-b4f2-7cf5c654d051
linearFit = GLM.predict(linearRegressor)

# ╔═╡ eeece605-6fb9-4a50-abd1-e7ab27602516
m = linearRegressor.model.pp.beta0[2]

# ╔═╡ a2ec1e1d-0cee-4239-8d44-e0072388da38
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 498ec62f-e555-46ce-b358-e1bb7a636c32
md"""
We can see from above that the linear regression returns a best fit line with the formula:

```math
y = mx + b
```

Which can be solved for ``x`` and then used to calculate the density (``x``) given some measured intensity (``y``)

```math
x = \frac{y - b}{m}
```

"""

# ╔═╡ e9c69ef8-7f37-483d-9704-f934a609f419
density(intensity) = (intensity - b) / m

# ╔═╡ 06136b05-d8b0-4cfe-8ec9-b6a64eff379a
intensity(ρ) = m * ρ + b

# ╔═╡ e09b2591-c9a8-4652-904a-9b2b7591afce
let
    f = Figure()
    ax1 = Axis(f[1, 1])

    scatter!(calibration_densities, calibration_intensities)
    lines!(calibration_densities, linearFit; color=:red)
    ax1.title = "Calibration Line (Intensity vs Density)"
    ax1.ylabel = "Intensity (HU)"
    ax1.xlabel = "Density (mg/cm^3)"

    f
end

# ╔═╡ 7b51613c-aa51-42af-8411-1fa37a4f909c
md"""
## Background mask
"""

# ╔═╡ 3c4b6370-14a4-4f95-9a19-22141666949b
ring_mask_L_HD = Bool.(dilate(dilate(dilate(dilate(mask_L_HD_3D)))) - dilate(dilate(dilate(mask_L_HD_3D))));

# ╔═╡ 39ad6549-867e-4161-ade6-747191e43c45
size(ring_mask_L_HD), size(dcm_array_inserts)

# ╔═╡ 73420cbd-9bf6-47d4-8528-1663cd2fc071
@bind g4 overlay_mask_bind(ring_mask_L_HD)

# ╔═╡ d3068023-5049-40c8-8fe6-339d44e65610
overlay_mask_plot(permutedims(dcm_array_inserts, (2, 1, 3)), permutedims(ring_mask_L_HD, (2, 1, 3)), g4, "Ring mask")

# ╔═╡ aef48319-82d4-483a-90b0-16efa7674d32
md"""
## Results
Since we have a calibration line, we don't actually need to directly measure the pure calcium signal ``S_{Obj}``. This is beneficial because segmenting pure calcium would be impractical for small calcifications. This also means we don't need to worry about heterogenous calcifications. We can assume the calcium is any density we want (within the calibration range), and if we keep the density consistent with the density used for ``S_{Obj}``, then the calculation will work out.
"""

# ╔═╡ b2949d94-c197-4c12-a96c-7e4dc736896d
S_Bkg = mean(dcm_array_inserts[ring_mask_L_HD])

# ╔═╡ 3f5a560c-0e03-4ca7-bb82-ba724cdb6742
S_Obj = intensity(400)

# ╔═╡ dfbc2297-3312-4ff5-85bc-6024c5447153
begin
    alg_integrated = Integrated(dcm_array_inserts[mask_L_HD_3D])
    ρ_hd = 0.4 # mg/mm^3
    mass_integrated_score = score(S_Bkg, S_Obj, spacing, ρ_hd, alg_integrated)
end

# ╔═╡ 42544384-6ff0-47ff-b071-d794cdf38022
md"""
We see that the ground truth mass (`mass_ground_truth` = $(mass_ground_truth)) is very similar to the calculated mass for integrated calcium scoring (`mass_integrated_score` = $(mass_integrated_score))
"""

# ╔═╡ 4ff23717-22fb-4b1b-998b-46b5c18996cf
md"""
# Spatially Weighted Calcium Scoring
Spatially weighted calcium scoring is another approach that attempts to improve the traditional Agatston approach. One way of doing this is by avoiding thresholding in favor of weighting each voxel. 

Below we will see how this technique can be used in CalciumScoring.jl
"""

# ╔═╡ e108fbb1-479a-4a2b-b279-abe502885df2
md"""
## Distribution
First, a distribution needs to be prepared based on measurements of 100 mg/cc calibration rods. This distribution can be estimated for the sake of documentation, but this should be measured carefully in practice.
"""

# ╔═╡ af011d0a-2c2f-429e-bbdd-595d2bd5ed14
μ, σ = 160, 30

# ╔═╡ 2d22a220-5849-439c-8705-234af4c219e5
md"""
## Results
Spatially weighted calcium scoring produces a score similar to Agatston scoring for high density calcifications, but more sensitive than Agatston scoring for low density calcifications. One limitation is that spatially weighted calcium scoring doesn't provide a straightforward way to convert from a spatially weighted calcium score to a physical measurement like volume score or mass score.
"""

# ╔═╡ 8e3c99eb-6be6-468e-a3ad-d2b6a69ae0d3
swcs_l_hd = score(overlayed_mask_l_hd, μ, σ, SpatiallyWeighted())

# ╔═╡ bb7f869d-0be7-4037-95b8-2d76ec540527
md"""
We see that the spatially weighted calcium score (`swcs_l_hd` = $(swcs_l_hd)) is very similar to the Agatston score (`agat_score` = $(agat_score))
"""

# ╔═╡ Cell order:
# ╟─c0bf04c7-7655-469a-8aaf-68f278bd75f6
# ╟─cc4ffe96-d49e-4ed2-955a-14b7d0617757
# ╠═2d99f63e-5489-11ed-33b9-4d6021e041dd
# ╠═ad36247d-9f0c-4083-bdee-21409b647c97
# ╟─830b0add-f62e-4340-8ee9-112f13e99326
# ╠═7c0763ec-cb65-40ae-9acf-339995c54284
# ╠═9525fdf5-cff1-41d2-b7b3-52f34df6af5a
# ╠═e4c99f57-f128-4759-bf58-02f5845d7d11
# ╠═45a05f47-f229-4d97-b5ae-929a27ac589e
# ╠═e1627af9-6dd7-4713-a9d1-e94cc7c4a890
# ╟─decb80cf-e696-47fd-9d50-04e0579b7fe6
# ╠═f4b2bd5f-dc8d-45ff-9629-a712c489ee60
# ╠═937bbe02-09c7-4c1e-8770-d5c681784958
# ╟─0afc60d7-9608-4189-b85a-7dc97a498e98
# ╟─6e6704d8-db30-4a17-ad08-c1928f5c1750
# ╟─83773d29-715f-45d4-a224-57d227c63f4a
# ╟─8cb8f78b-515c-46ee-9edc-ea22b538a5b5
# ╟─bb768658-c623-4941-b248-30fb6fe76779
# ╟─6ad2569c-2e98-4abe-a02e-80943e8b7b2f
# ╟─b8738edd-3be2-44fc-9812-a96332c51f5f
# ╠═f07faa8a-6653-46ac-93b1-8e8d6ef53a94
# ╠═89b45727-74ab-4c3b-bf24-53642def397f
# ╠═fca9ac32-e940-4a6a-8f84-51547c5329ce
# ╠═e27ed808-20cf-4846-b2d3-53ef3dd9ff71
# ╟─d453b50f-06e5-4cab-a6b7-a0f1629770a9
# ╟─6bad9fee-b5b1-4c58-a77c-95cff35982ef
# ╠═480dcaca-861a-4d51-8bcc-dc3e24980e9c
# ╟─1fafecbf-b767-49dc-938f-ad2b6c4241de
# ╠═555eac6e-ce5c-4002-8a31-7e39ed0174c3
# ╟─00fcc622-cf8b-43ee-a5bd-acafabc903f7
# ╟─a4510ff6-92df-4aac-a774-7f046cfa9ea4
# ╟─dff1746c-6e43-4cdc-a41d-b3c26ee9e680
# ╠═a2d5117a-04e2-4476-9ef6-6d38070a2264
# ╠═682c3d99-11c6-4927-bcf5-0c1833a5588e
# ╠═b26f3244-ea93-41aa-b4f2-7cf5c654d051
# ╠═eeece605-6fb9-4a50-abd1-e7ab27602516
# ╠═a2ec1e1d-0cee-4239-8d44-e0072388da38
# ╟─498ec62f-e555-46ce-b358-e1bb7a636c32
# ╠═e9c69ef8-7f37-483d-9704-f934a609f419
# ╠═06136b05-d8b0-4cfe-8ec9-b6a64eff379a
# ╟─e09b2591-c9a8-4652-904a-9b2b7591afce
# ╟─7b51613c-aa51-42af-8411-1fa37a4f909c
# ╠═3c4b6370-14a4-4f95-9a19-22141666949b
# ╠═39ad6549-867e-4161-ade6-747191e43c45
# ╟─73420cbd-9bf6-47d4-8528-1663cd2fc071
# ╟─d3068023-5049-40c8-8fe6-339d44e65610
# ╟─aef48319-82d4-483a-90b0-16efa7674d32
# ╠═b2949d94-c197-4c12-a96c-7e4dc736896d
# ╠═3f5a560c-0e03-4ca7-bb82-ba724cdb6742
# ╠═dfbc2297-3312-4ff5-85bc-6024c5447153
# ╟─42544384-6ff0-47ff-b071-d794cdf38022
# ╟─4ff23717-22fb-4b1b-998b-46b5c18996cf
# ╟─e108fbb1-479a-4a2b-b279-abe502885df2
# ╠═af011d0a-2c2f-429e-bbdd-595d2bd5ed14
# ╟─2d22a220-5849-439c-8705-234af4c219e5
# ╠═8e3c99eb-6be6-468e-a3ad-d2b6a69ae0d3
# ╟─bb7f869d-0be7-4037-95b8-2d76ec540527

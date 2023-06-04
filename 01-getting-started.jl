### A Pluto.jl notebook ###
# v0.19.26

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

# ╔═╡ 9419f973-985f-41af-b6a4-bbd6af9ea2a4
# ╠═╡ show_logs = false
begin
	let
		using Pkg
		Pkg.activate(temp = true)
		Pkg.Registry.update()
		Pkg.add("PlutoUI")
		Pkg.add("CairoMakie")
		Pkg.add("CSV")
		Pkg.add("DataFrames")
		Pkg.add("ImageMorphology")
		Pkg.add("GLM")
		Pkg.add("Statistics")
		Pkg.add("ImageFiltering")
		Pkg.add("Noise")
		Pkg.add("DICOM")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI, CSV, DataFrames, CairoMakie, ImageMorphology, GLM, Statistics, ImageFiltering, Noise, DICOM, CalciumScoring
end

# ╔═╡ d6254d97-3b3b-4b2f-b7b3-8bc58e6f34c2
md"""
# Overview
This tutorial will briefly showcase different calcium scoring techniques available in CalciumScoring.jl. Specifically (1) **Agatston Scoring** which is the current gold standard, (2) **integrated calcium mass**, (3) **volume fraction calcium mass**, and (4) **spatially weighted calcium scoring**

These techniques provide quantitative measures of calcium from non-contrast CT images. Every calcium mass technique calculates calcium mass via a calibration factor.

Calcium scoring also contains **material decomposition calcium mass**, another calcium measurement technique based on dual-energy material decomposition. [Click here](https://glassnotebook.io/repo/31) to get started using the material decomposition calcium mass technique.
"""

# ╔═╡ f3a886e7-31e8-43dd-adc7-3fb50c7d7d23
md"""
## Import Packages
First, let's import the most up-to-date version of CalciumScoring.jl, which can be found on the main/master branch of the [GitHub repository](https://github.com/Dale-Black/CalciumScoring.jl).

Because we are using the unregistered version (most recent), we will need to `Pkg.add` this explicitly without using Pluto's built-in package manager. Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.

To help with the formatting of this documentation, we will also add [PlutoUI.jl](https://github.com/JuliaPluto/PlutoUI.jl). We will also create some simulated coronary artery images, which will require [ImageFiltering.jl](https://github.com/JuliaImages/ImageFiltering.jl) and [Noise.jl](https://github.com/roflmaostc/Noise.jl) Lastly, to visualize the results, we're going to add [Makie.jl](https://github.com/JuliaPlots/Makie.jl).
"""

# ╔═╡ df5551ce-01cc-49a2-b040-e1a2b6d794a4
TableOfContents()

# ╔═╡ 87286995-8337-410d-bc12-a2c1c638aff5
md"""
## Calculate Ground Truth Values
Since we know the intensity value for every calcium voxel in the `pure_calcification_3D` array should be 200, we can calculate the number of calcium voxels in that array directly. We can then multiply this by the voxel size to determine the ground truth calcium volume. Lastly, we can take the density of calcium and multiply it by the volume to determine the mass of the calcium
"""

# ╔═╡ b0afd734-56ea-45df-bf68-a282a9ef7112
ground_truth_number_voxels = length(findall(x -> x == hu_calcium, pure_calcification_3D))

# ╔═╡ b667c99e-a012-4418-b641-b7dfce71332f
ground_truth_calcium_volume = ground_truth_number_voxels * voxel_size

# ╔═╡ 93123baf-d4ad-4959-b4df-7f859654f191
ground_truth_calcium_mass = ground_truth_calcium_volume * ρ_calcium

# ╔═╡ 0188bed1-93f4-4e95-b581-0ea64847faac
md"""
# Visualize Calcification
Below, in the appendix, we created a simulated image containing a calcification at the center. Let's visualize this image below.
"""

# ╔═╡ 5c64de6c-2143-43a5-94cf-a7cdb6e4e210
@bind x PlutoUI.Slider(axes(gaussian_poisson_calcification_3D, 3); default=2, show_value=true)

# ╔═╡ 0a3c60e5-b695-41e4-b43d-45b2b56c5b7e
let
	f = Figure()
	ax = CairoMakie.Axis(
		f[1, 1],
		title = "Simulated Calcification Image"
	)
	heatmap!(gaussian_poisson_calcification_3D[:, :, x], colormap=:grays)
	f
end

# ╔═╡ 13bf4d31-e155-41ee-81c1-55b2b6f28ecd
md"""
# 1. Agatston Scoring
Calcium scoring is a technique for measuring calcium in the coronary arteries within a CT scan. Calcium scoring is traditionally associated with the Agatston scoring method [(1)](https://pubmed.ncbi.nlm.nih.gov/2407762/)

This package allows users to compute the traditional Agatston score seen below.
"""

# ╔═╡ 31ba16d9-e51f-48ef-8731-fce208de1c30
md"""
## Visualize
Now let's dilate the mask a little bit, to account for partial volume effect (blurring) and visualize our ROI
"""

# ╔═╡ 51649204-97fe-460b-935e-06a29bb85d5b
dilated_mask_3D = dilate(dilate(dilate(dilate(dilate(dilate(dilate(mask_3D)))))));

# ╔═╡ 943ad7b3-965a-45e9-b135-ad40432d1fde
md"""
## Results
"""

# ╔═╡ 8b13e40c-b2fd-4f2d-856e-7605c5f5c0b7
md"""
We can also compute the calcium mass score, via the Agatston technique by computing the mass calibration factor. This factor translates an Agatston score to a mass score ([1](https://pubs.rsna.org/doi/epdf/10.1148/radiol.2432050808)).

We can assume the CT number of water is 0, and compute the correction factor by dividing the density of known calcium by the CT number of that same calcium.
"""

# ╔═╡ c753ffb7-fef1-45db-b43a-e368b2fd1cdf
md"""
# 2. Integrated Calcium Mass
For integrated calcium mass, we need the dilated mask from above, along with a background mask, and a calibration line for ``S_{Obj}``
"""

# ╔═╡ 919d0be1-630a-41b1-9218-3a69496acffb
md"""
## Calibration Line
For a true calibration line, one would need to include known calcium density inserts in the original image for patient-specific calibration. For this example, we will use previously validated calibration points, specific to these simulated images. This will allow us to then compute a calibration line, using GLM.jl
"""

# ╔═╡ 477b8213-be8a-413b-a9cb-4bf5c3492747
begin
	calibration_densities = [0, 0.2, 0.4, 0.8]mg/mm^3
	calibration_intensities = [0, 297.429, 545.245, 1075.82] # HU
	df_calibration = DataFrame(:density => calibration_densities, :intensity => calibration_intensities)
end

# ╔═╡ 60b0f243-468c-4368-b9f0-d734edd476c0
linearRegressor = lm(@formula(intensity ~ density), ustrip.(df_calibration));

# ╔═╡ 2e35e337-e64c-4a62-af2d-bf93151781c2
linearFit = GLM.predict(linearRegressor)

# ╔═╡ a262682b-798a-4531-8281-a6acfb9560f4
m = linearRegressor.model.pp.beta0[2]mm^3/mg

# ╔═╡ 4889b128-4bfb-419d-99ac-4bb3ef0ebba6
b = linearRegressor.model.rr.mu[1]

# ╔═╡ 147d9712-d861-4182-a44c-2d84ebd82200
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

# ╔═╡ ac3e5c8e-6006-4b33-8711-24792c64a463
intensity(ρ) = m * ρ + b

# ╔═╡ cc2a8b81-8506-45bc-bdc9-c5f1aa7d6c83
let
    f = Figure()
    ax1 = Axis(f[1, 1])

    scatter!(ustrip.(calibration_densities), calibration_intensities)
    lines!(ustrip.(calibration_densities), linearFit; color=:red)
    ax1.title = "Calibration Line (Intensity vs Density)"
    ax1.ylabel = "Intensity (HU)"
    ax1.xlabel = "Density (mg/cm^3)"

    f
end

# ╔═╡ ae1cd5d6-ecee-44c6-aaf6-8bc314fd737c
md"""
## Background mask
"""

# ╔═╡ 0cb27e66-f2c0-4245-898c-9d2ebdcd0636
ring_mask_3D = Bool.(dilate(dilate(dilate(dilate(dilate(dilate(dilate(mask_3D))))))) - dilate(dilate(dilate(dilate(dilate(mask_3D))))));

# ╔═╡ 8fe84415-a1e2-4180-bc4d-02dc1b3fbd08
size(ring_mask_3D), size(gaussian_poisson_calcification_3D)

# ╔═╡ 205d0159-13ef-4483-98c5-b4f86bb59510
gaussian_poisson_calcification_3D[ring_mask_3D]

# ╔═╡ ffff8cbd-fbeb-4a5d-bb3f-560ec2973137
md"""
## Results
Since we have a calibration line, we don't actually need to directly measure the pure calcium signal ``S_{Obj}``. This is beneficial because segmenting pure calcium would be impractical for small calcifications. This also means we don't need to worry about heterogenous calcifications. We can assume the calcium is any density we want (within the calibration range), and if we keep the density consistent with the density used for ``S_{Obj}``, then the calculation will work out.
"""

# ╔═╡ f63e8885-aa22-418c-b3db-886fe1f2945e
S_Bkg = mean(gaussian_poisson_calcification_3D[ring_mask_3D])

# ╔═╡ e574fb88-fa85-4efc-b64b-c1bb39473cec
	S_Obj = intensity(ρ_calcium)

# ╔═╡ 04446f5f-eda2-4ae9-93e8-ccba84b26014
begin
    alg_integrated = Integrated(gaussian_poisson_calcification_3D[dilated_mask_3D])
    mass_integrated_score = score(S_Bkg, S_Obj, spacing, ρ_calcium, alg_integrated)
end

# ╔═╡ 128e32e5-4ddb-4f18-875f-e28e5b548d2c
md"""
Let's compare that with the ground truth mass

We see that the ground truth mass = $(ground_truth_calcium_mass) is close to the calculated mass = $(mass_integrated_score)
"""

# ╔═╡ 122317ec-7c00-4b1e-9045-1314f444e56d
md"""
# 3. Volume Fraction Calcium Mass
We will examine a simple yet powerful calcium quantification technique. This technique is similar to the integrated calcium mass technique, with fewer steps. All that is required is a region of interest containing all the calcium, the known HU of a specific calcium calibration rod, and the known HU of background material.
"""

# ╔═╡ 017fe43a-df10-40d9-b9d1-d93d71500bd0
md"""
## Results
"""

# ╔═╡ 7635fbef-61fa-4311-bda4-f7f384a4832d
volume_fraction_mass = score(gaussian_poisson_calcification_3D[dilated_mask_3D], hu_calcium, hu_heart_tissue, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 11c49d54-2088-4f6c-83a1-129960889492
md"""
Let's compare that with the ground truth mass

We see that the ground truth mass = $(ground_truth_calcium_mass) is close to the calculated mass = $(volume_fraction_mass)
"""

# ╔═╡ a91fd5a5-cd7b-4676-a6c3-13d8c81453d9
md"""
# 4. Spatially Weighted Calcium Scoring
Spatially weighted calcium scoring is another approach that attempts to improve upon the traditional Agatston approach. One way of doing this is by avoiding thresholding in favor of weighting each voxel. 

Below we will see how this technique can be used in CalciumScoring.jl
"""

# ╔═╡ c624e494-ade4-4a05-98e7-f7beaeaf847f
md"""
## Distribution
First, a distribution needs to be prepared based on measurements of 100 mg/cc calibration rods. This distribution can be estimated for the sake of documentation, but this should be measured carefully in practice.
"""

# ╔═╡ a0609e90-45ea-4494-ad0a-86d339dff917
μ, σ = 175, 20

# ╔═╡ 581e9601-c169-4b2a-acec-391a75be26da
md"""
## Results
Spatially weighted calcium scoring produces a score similar to Agatston scoring for high density calcifications, but more sensitive than Agatston scoring for low density calcifications. One limitation is that spatially weighted calcium scoring doesn't provide a straightforward way to convert from a spatially weighted calcium score to a physical measurement like volume score or mass score.
"""

# ╔═╡ badbb0fb-6034-4111-87c6-a0188b73ea14
md"""
# Results Summary
"""

# ╔═╡ 7484f5a6-0cbc-451f-b790-aea05813c425
md"""
# Appendix
"""

# ╔═╡ fac54db6-e2a1-4956-8b7c-091938e74b8e
md"""
#### Load Simulated Phantoms
"""

# ╔═╡ 29b04989-caa7-4e73-a9a1-cdad9e3509d8
begin
	dcm_path = joinpath(pwd(), "dcms/val/52_59_73/medium/80");
	dcms = dcmdir_parse(dcm_path)
	dcm_array = cat([dcms[i][(0x7fe0, 0x0010)] for i in 1:length(dcms)]...; dims=3)
end;

# ╔═╡ 982837e3-6d4c-4012-a014-886d38d19121
md"""
#### Helper functions
"""

# ╔═╡ 107e9699-6baa-4978-a68f-9fc22d76d6a1
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ ac9da8e2-aee9-4c79-8e7d-f35eb7e674f2
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=2, show_value=true)
end

# ╔═╡ 83246419-b0e2-478c-bb0d-a7ab12c6ecd7
@bind c overlay_mask_bind(dilated_mask_3D)

# ╔═╡ f0b9b969-eca7-4da8-a95f-a5aef1597a32
@bind g4 overlay_mask_bind(ring_mask_3D)

# ╔═╡ 4a14052c-5327-43a4-bf2e-ad739bf5213d
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
        markersize=1.2,
        color=:red,
    )
    return fig
end

# ╔═╡ 907fe9d2-a67e-4716-9c73-0ce429ae2b4f
overlay_mask_plot(gaussian_poisson_calcification_3D, dilated_mask_3D, c, "Overlayed Mask")

# ╔═╡ 98448e7c-3c24-4ede-9adb-3691e68d0742
overlay_mask_plot(gaussian_poisson_calcification_3D, ring_mask_3D, g4, "Ring mask")

# ╔═╡ 140d7064-df9e-4dc8-886e-edb381792164
function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

# ╔═╡ e77c3179-ade5-47d2-b81c-c97087161632
overlayed_mask = create_mask(gaussian_poisson_calcification_3D, dilated_mask_3D);

# ╔═╡ 913619e8-80ae-4130-831e-9f12941e16a2
heatmap(overlayed_mask[:, :, c], colormap=:grays)

# ╔═╡ ade72d79-4082-4c4e-b14a-744651731d9d
begin
	mass_calibration = ρ_calcium / hu_calcium
	agatston_score, volume_score, mass_score = score(overlayed_mask, spacing, mass_calibration, Agatston())
end

# ╔═╡ 320405f5-4feb-4c61-98f4-ff13643490ae
md"""
Let's compare that with the ground truth mass

We see that the ground truth mass = $(ground_truth_calcium_mass) is close to the calculated mass = $(mass_score)
"""

# ╔═╡ 8df54213-6240-45d5-baaf-36a81094f998
calculated_masses = [mass_score, mass_integrated_score, volume_fraction_mass]

# ╔═╡ 879aafa3-a9b3-4a8e-bd54-6068be270332
df_results_masses = DataFrame(
	technique = ["Agatston", "Integrated", "Volume Fraction"],
	ground_truth_calcium_mass = repeat([ground_truth_calcium_mass], 3),
	calculated_calcium_mass = calculated_masses
)

# ╔═╡ 15a6ff92-a3c1-4ff1-bc4b-674adbc53c0d
swcs = score(overlayed_mask, μ, σ, SpatiallyWeighted())

# ╔═╡ 1e763d9c-44e8-4614-94fb-be9b65f100a9
md"""
We see that the spatially weighted calcium score (`swcs` = $(swcs)) is very similar to the Agatston score (`agatston_score` = $(agatston_score))
"""

# ╔═╡ 757901be-3d2f-4bea-9054-469cf7a18998
calculated_scores = [agatston_score, swcs]

# ╔═╡ f46b4cd0-cfc1-48fd-b964-b8cc282f8f53
df_results = DataFrame(
	technique = ["Agatston", "Spatially Weighted"],
	calculated_scores = calculated_scores
)

# ╔═╡ Cell order:
# ╟─d6254d97-3b3b-4b2f-b7b3-8bc58e6f34c2
# ╟─f3a886e7-31e8-43dd-adc7-3fb50c7d7d23
# ╠═9419f973-985f-41af-b6a4-bbd6af9ea2a4
# ╠═df5551ce-01cc-49a2-b040-e1a2b6d794a4
# ╟─87286995-8337-410d-bc12-a2c1c638aff5
# ╠═b0afd734-56ea-45df-bf68-a282a9ef7112
# ╠═b667c99e-a012-4418-b641-b7dfce71332f
# ╠═93123baf-d4ad-4959-b4df-7f859654f191
# ╟─0188bed1-93f4-4e95-b581-0ea64847faac
# ╟─5c64de6c-2143-43a5-94cf-a7cdb6e4e210
# ╟─0a3c60e5-b695-41e4-b43d-45b2b56c5b7e
# ╟─13bf4d31-e155-41ee-81c1-55b2b6f28ecd
# ╟─31ba16d9-e51f-48ef-8731-fce208de1c30
# ╠═51649204-97fe-460b-935e-06a29bb85d5b
# ╟─83246419-b0e2-478c-bb0d-a7ab12c6ecd7
# ╟─907fe9d2-a67e-4716-9c73-0ce429ae2b4f
# ╠═e77c3179-ade5-47d2-b81c-c97087161632
# ╟─913619e8-80ae-4130-831e-9f12941e16a2
# ╟─943ad7b3-965a-45e9-b135-ad40432d1fde
# ╟─8b13e40c-b2fd-4f2d-856e-7605c5f5c0b7
# ╠═ade72d79-4082-4c4e-b14a-744651731d9d
# ╟─320405f5-4feb-4c61-98f4-ff13643490ae
# ╟─c753ffb7-fef1-45db-b43a-e368b2fd1cdf
# ╟─919d0be1-630a-41b1-9218-3a69496acffb
# ╠═477b8213-be8a-413b-a9cb-4bf5c3492747
# ╠═60b0f243-468c-4368-b9f0-d734edd476c0
# ╠═2e35e337-e64c-4a62-af2d-bf93151781c2
# ╠═a262682b-798a-4531-8281-a6acfb9560f4
# ╠═4889b128-4bfb-419d-99ac-4bb3ef0ebba6
# ╟─147d9712-d861-4182-a44c-2d84ebd82200
# ╠═ac3e5c8e-6006-4b33-8711-24792c64a463
# ╟─cc2a8b81-8506-45bc-bdc9-c5f1aa7d6c83
# ╟─ae1cd5d6-ecee-44c6-aaf6-8bc314fd737c
# ╠═0cb27e66-f2c0-4245-898c-9d2ebdcd0636
# ╠═8fe84415-a1e2-4180-bc4d-02dc1b3fbd08
# ╟─f0b9b969-eca7-4da8-a95f-a5aef1597a32
# ╠═98448e7c-3c24-4ede-9adb-3691e68d0742
# ╠═205d0159-13ef-4483-98c5-b4f86bb59510
# ╟─ffff8cbd-fbeb-4a5d-bb3f-560ec2973137
# ╠═f63e8885-aa22-418c-b3db-886fe1f2945e
# ╠═e574fb88-fa85-4efc-b64b-c1bb39473cec
# ╠═04446f5f-eda2-4ae9-93e8-ccba84b26014
# ╟─128e32e5-4ddb-4f18-875f-e28e5b548d2c
# ╟─122317ec-7c00-4b1e-9045-1314f444e56d
# ╟─017fe43a-df10-40d9-b9d1-d93d71500bd0
# ╠═7635fbef-61fa-4311-bda4-f7f384a4832d
# ╟─11c49d54-2088-4f6c-83a1-129960889492
# ╟─a91fd5a5-cd7b-4676-a6c3-13d8c81453d9
# ╟─c624e494-ade4-4a05-98e7-f7beaeaf847f
# ╠═a0609e90-45ea-4494-ad0a-86d339dff917
# ╟─581e9601-c169-4b2a-acec-391a75be26da
# ╠═15a6ff92-a3c1-4ff1-bc4b-674adbc53c0d
# ╟─1e763d9c-44e8-4614-94fb-be9b65f100a9
# ╟─badbb0fb-6034-4111-87c6-a0188b73ea14
# ╠═8df54213-6240-45d5-baaf-36a81094f998
# ╠═879aafa3-a9b3-4a8e-bd54-6068be270332
# ╠═757901be-3d2f-4bea-9054-469cf7a18998
# ╠═f46b4cd0-cfc1-48fd-b964-b8cc282f8f53
# ╟─7484f5a6-0cbc-451f-b790-aea05813c425
# ╟─fac54db6-e2a1-4956-8b7c-091938e74b8e
# ╠═29b04989-caa7-4e73-a9a1-cdad9e3509d8
# ╟─982837e3-6d4c-4012-a014-886d38d19121
# ╟─107e9699-6baa-4978-a68f-9fc22d76d6a1
# ╟─ac9da8e2-aee9-4c79-8e7d-f35eb7e674f2
# ╟─4a14052c-5327-43a4-bf2e-ad739bf5213d
# ╟─140d7064-df9e-4dc8-886e-edb381792164

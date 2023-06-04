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
This tutorial will briefly showcase different calcium scoring techniques available in CalciumScoring.jl. Specifically (1) **Agatston Scoring** which is the current gold standard, (2) **volume fraction calcium mass**, and (3) **spatially weighted calcium scoring**

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

# ╔═╡ 0188bed1-93f4-4e95-b581-0ea64847faac
md"""
# Visualize Simulated Phantom
Below, in the appendix, we created a simulated phantom containing a calcification at the center. Let's visualize this image below.

This phantom is supposed to simulate a human thorax containing calcified coronary arteries. Since the phantom has identical geometry along each slice, the only difference along each slice will be from "quantum noise".
"""

# ╔═╡ 87286995-8337-410d-bc12-a2c1c638aff5
md"""
## Calculate Ground Truth Values
The goal with calcium scoring is to quantify the amount of calcium within a cardiac CT scan. Before we calculate the calcium via various calcium scoring techniques, lets calculate the know calcium mass. Since we know the density of the calcifications and the geometry of the inserts, we can calculate the ground truth calcium mass of each of the 9 calcium inserts. 

The calculation for volume is simple since the inserts are cylindrical
```math
\text{volume} = \pi \times r^2 \times \text{slice thickness} \times \text{number of slices}
\tag{1}
```

Mass is then simple to obtain from the volume and known density
```math
\text{mass} = \text{volume} \times \rho
\tag{1}
```

"""

# ╔═╡ ae4afb28-3caf-4615-8e89-774dcef8dd7e
ρ_high_density = 73 # mgcm-3

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

# ╔═╡ 310ad8ab-838a-4080-ad09-076a56767174
md"""
## Results
"""

# ╔═╡ 8b13e40c-b2fd-4f2d-856e-7605c5f5c0b7
md"""
We can also compute the calcium mass score, via the Agatston technique by computing the mass calibration factor. This factor translates an Agatston score to a mass score ([1](https://pubs.rsna.org/doi/epdf/10.1148/radiol.2432050808)).

We can assume the CT number of water is 0, and compute the correction factor by dividing the density of known calcium by the CT number of that same calcium. Let's load our simulated calibration phantom to acquire this.
"""

# ╔═╡ c3837d06-7c00-49ca-8592-6b9e8e7ef528
center_insert1, center_insert2 = 187, 318

# ╔═╡ 122317ec-7c00-4b1e-9045-1314f444e56d
md"""
# 2. Volume Fraction Calcium Mass
We will examine a simple yet powerful calcium quantification technique. This technique is similar to the integrated calcium mass technique, with fewer steps. All that is required is a region of interest containing all the calcium, the known HU of a specific calcium calibration rod, and the known HU of background material.
"""

# ╔═╡ 017fe43a-df10-40d9-b9d1-d93d71500bd0
md"""
## Results
"""

# ╔═╡ a91fd5a5-cd7b-4676-a6c3-13d8c81453d9
md"""
# 3. Spatially Weighted Calcium Scoring
Spatially weighted calcium scoring is another approach that attempts to improve upon the traditional Agatston approach. One way of doing this is by avoiding thresholding in favor of weighting each voxel. 

Below we will see how this technique can be used in CalciumScoring.jl
"""

# ╔═╡ c624e494-ade4-4a05-98e7-f7beaeaf847f
md"""
#### Distribution
First, a distribution needs to be prepared based on measurements of 100 mg/cc calibration rods. This distribution can be estimated for the sake of documentation, but this should be measured carefully in practice.
"""

# ╔═╡ a0609e90-45ea-4494-ad0a-86d339dff917
μ, σ = 175, 20

# ╔═╡ 581e9601-c169-4b2a-acec-391a75be26da
md"""
## Results
Spatially weighted calcium scoring produces a score similar to Agatston scoring for high density calcifications, but more sensitive than Agatston scoring for low density calcifications. One limitation is that spatially weighted calcium scoring doesn't provide a straightforward way to convert from a spatially weighted calcium score to a physical measurement like volume score or mass score.
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
	# Validation Phantom
	dcm_path = joinpath(pwd(), "dcms/val/52_59_73/medium/80");
	dcms = dcmdir_parse(dcm_path)
	dcm_array = cat([dcms[i][(0x7fe0, 0x0010)] for i in 1:length(dcms)]...; dims=3)

	# Calibration Phantom
	dcm_path_calibration = joinpath(pwd(), "dcms/cal/200/30/80");
	dcms_calibration = dcmdir_parse(dcm_path_calibration)
	dcm_array_calibration = cat([dcms_calibration[i][(0x7fe0, 0x0010)] for i in 1:length(dcms_calibration)]...; dims=3)

	# Calcification Masks
	masks_medium = Dict(
		:mask_L_HD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks", "medium", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_M_HD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_S_HD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_L_MD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_M_MD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_S_MD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_M_LD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_L_LD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_S_LD => Array(CSV.read(joinpath(pwd(), "dcms", "val_masks","medium", "mask_S_LD.csv"), DataFrame; header=false)),
	)
end;

# ╔═╡ 5c64de6c-2143-43a5-94cf-a7cdb6e4e210
@bind z1 PlutoUI.Slider(axes(dcm_array, 3); default = 2, show_value = true)

# ╔═╡ 9c700f7b-53b9-4f03-bef1-fef6a0e3cc68
slice_thickness = dcms[1].meta[tag"Slice Thickness"]

# ╔═╡ b25cf073-a3ab-4704-86ca-18f1ecbf6418
number_of_slices = length(axes(dcm_array, 3))

# ╔═╡ b667c99e-a012-4418-b641-b7dfce71332f
begin
	ground_truth_calcium_volume_large_mm3 = π * (5/2)^2 * slice_thickness * number_of_slices # mm^3
	ground_truth_calcium_volume_large_cm3 = ground_truth_calcium_volume_large_mm3 * 1e-3
end

# ╔═╡ 93123baf-d4ad-4959-b4df-7f859654f191
ground_truth_calcium_mass_mg_large_high_density = ground_truth_calcium_volume_large_cm3 * ρ_high_density

# ╔═╡ 83246419-b0e2-478c-bb0d-a7ab12c6ecd7
@bind z2 PlutoUI.Slider(axes(dcm_array, 3); show_value = true, default = 2)

# ╔═╡ 4d8089c3-435e-4b1a-a185-bd79e1425c13
pixel_spacing = dcms[1][tag"Pixel Spacing"]

# ╔═╡ 9fc7b12c-b770-4431-bc38-0293b88a13f1
voxel_spacing = [pixel_spacing..., slice_thickness]

# ╔═╡ d37b1618-f0ed-42af-a336-8f9a341f2445
voxel_size = voxel_spacing[1] * voxel_spacing[2] * voxel_spacing[1]

# ╔═╡ f00fb7d4-5183-483a-90a2-73ff9be3a24b
@bind z3 PlutoUI.Slider(axes(dcm_array_calibration, 3); default = 2, show_value = true)

# ╔═╡ efd98ea7-b86a-41e1-9eca-6675dffc7edd
let
	
	f = Figure()
	ax = Axis(
		f[1, 1],
		title = "DICOM Calibration Rod"
	)
	heatmap!(transpose(dcm_array_calibration[:, :, z2]); colormap = :grays)
	f
end

# ╔═╡ 1870556b-2e19-4e64-9df1-338daae29805
begin
	calibration_rod = zeros(25, 25, size(dcm_array_calibration, 3))
				
	for z in axes(dcm_array_calibration, 3)
		rows, cols, depth = size(dcm_array_calibration)
		half_row, half_col = center_insert1, center_insert2
		offset = 12
		row_range = half_row-offset:half_row+offset
		col_range = half_col-offset:half_col+offset	
		calibration_rod[:, :, z] .= dcm_array_calibration[row_range, col_range, z];
	end
	hu_calcium_200, ρ_calcium_200 = mean(calibration_rod), 200
end

# ╔═╡ 2d3b75aa-da53-4490-a743-9a6ca405759f
mass_calibration = ρ_calcium_200 / hu_calcium_200 * 1e-3

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

# ╔═╡ 78d99536-b253-491e-b08d-df7cb97a2063
"""
    dilate_recursively(mask, n)

Recursively erode a `mask`

#### Inputs
- `mask`: boolean array to erode
- `n`: number of erosions

#### Returns
- `dilated_mask`: dilated mask with `n` erosions
"""
function dilate_recursively(mask, n)
    dilated_mask = copy(mask)
    for _ in 1:n
        dilated_mask = dilate(dilated_mask)
    end
    return dilated_mask
end;

# ╔═╡ 51649204-97fe-460b-935e-06a29bb85d5b
begin
	dilated_mask_large_high_density = dilate_recursively(masks_medium[:mask_L_LD], 2)
	dilated_mask_large_high_density_3D = cat(
		dilated_mask_large_high_density, 
		dilated_mask_large_high_density,
		dilated_mask_large_high_density;
		dims=3)
end;

# ╔═╡ 0a3c60e5-b695-41e4-b43d-45b2b56c5b7e
let
	idxs = Tuple.(findall(isone, dilated_mask_large_high_density_3D[:, :, z1]))
	idxs = [t[i] for t in idxs, i in 1:length(idxs[1])]
	
	f = Figure(resolution = (1200, 800))
	ax = Axis(
		f[1, 1],
		title = "DICOM Array"
	)
	heatmap!(transpose(dcm_array[:, :, z1]); colormap = :grays)
	f
end

# ╔═╡ 41340d94-d8e6-48cd-879f-707017969911
pure_calcification_large_high_density = create_mask(dcm_array, dilated_mask_large_high_density_3D);

# ╔═╡ 7e4c3269-2a2b-4b15-a373-4daf815a72ef
agatston_score, volume_score = score(pure_calcification_large_high_density, voxel_spacing, Agatston())

# ╔═╡ ade72d79-4082-4c4e-b14a-744651731d9d
_, _, mass_score = score(pure_calcification_large_high_density, voxel_spacing, mass_calibration, Agatston())

# ╔═╡ 320405f5-4feb-4c61-98f4-ff13643490ae
md"""
Let's compare that with the ground truth mass

We see that the ground truth mass = $(ground_truth_calcium_mass_mg_large_high_density) is close to the calculated mass = $(mass_score)
"""

# ╔═╡ 15a6ff92-a3c1-4ff1-bc4b-674adbc53c0d
swcs = score(pure_calcification_large_high_density, μ, σ, SpatiallyWeighted())

# ╔═╡ 1e763d9c-44e8-4614-94fb-be9b65f100a9
md"""
We see that the spatially weighted calcium score (`swcs` = $(swcs)) is very similar to the Agatston score (`agatston_score` = $(agatston_score))
"""

# ╔═╡ 907fe9d2-a67e-4716-9c73-0ce429ae2b4f
let
	idxs = Tuple.(findall(isone, dilated_mask_large_high_density_3D[:, :, z2]))
	idxs = [t[i] for t in idxs, i in 1:length(idxs[1])]
	
	f = Figure(resolution = (1200, 800))
	
	ax = Axis(
		f[1, 1],
		title = "Overlayed Mask"
	)
	heatmap!(transpose(dcm_array[:, :, z2]); colormap = :grays)
	scatter!(idxs[:, 2], idxs[:, 1]; markersize = 1, color = :red)

	ax = Axis(
		f[1, 2],
		title = "Pure Calcificaiton"
	)
	heatmap!(transpose(pure_calcification_large_high_density[:, :, z2]))
	f
end

# ╔═╡ 7635fbef-61fa-4311-bda4-f7f384a4832d
volume_fraction_mass = score(dcm_array[dilated_mask_large_high_density_3D], hu_calcium_200, 30, voxel_size, ρ_calcium_200, VolumeFraction())

# ╔═╡ 11c49d54-2088-4f6c-83a1-129960889492
md"""
Let's compare that with the ground truth mass

We see that the ground truth mass = $(ground_truth_calcium_mass_mg_large_high_density) is close to the calculated mass = $(volume_fraction_mass)
"""

# ╔═╡ 7778b027-68be-46c7-a8f3-31413fa5554f
"""
    erode_recursively(mask, n)

Recursively erode a `mask`

#### Inputs
- `mask`: boolean array to erode
- `n`: number of erosions

#### Returns
- `eroded_mask`: eroded mask with `n` erosions
"""
function erode_recursively(mask, n)
    eroded_mask = copy(mask)
    for _ in 1:n
        eroded_mask = erode(eroded_mask)
    end
    return eroded_mask
end;

# ╔═╡ Cell order:
# ╟─d6254d97-3b3b-4b2f-b7b3-8bc58e6f34c2
# ╟─f3a886e7-31e8-43dd-adc7-3fb50c7d7d23
# ╠═9419f973-985f-41af-b6a4-bbd6af9ea2a4
# ╠═df5551ce-01cc-49a2-b040-e1a2b6d794a4
# ╟─0188bed1-93f4-4e95-b581-0ea64847faac
# ╟─5c64de6c-2143-43a5-94cf-a7cdb6e4e210
# ╟─0a3c60e5-b695-41e4-b43d-45b2b56c5b7e
# ╟─87286995-8337-410d-bc12-a2c1c638aff5
# ╠═9c700f7b-53b9-4f03-bef1-fef6a0e3cc68
# ╠═b25cf073-a3ab-4704-86ca-18f1ecbf6418
# ╠═ae4afb28-3caf-4615-8e89-774dcef8dd7e
# ╠═b667c99e-a012-4418-b641-b7dfce71332f
# ╠═93123baf-d4ad-4959-b4df-7f859654f191
# ╟─13bf4d31-e155-41ee-81c1-55b2b6f28ecd
# ╟─31ba16d9-e51f-48ef-8731-fce208de1c30
# ╠═51649204-97fe-460b-935e-06a29bb85d5b
# ╠═41340d94-d8e6-48cd-879f-707017969911
# ╟─83246419-b0e2-478c-bb0d-a7ab12c6ecd7
# ╟─907fe9d2-a67e-4716-9c73-0ce429ae2b4f
# ╟─310ad8ab-838a-4080-ad09-076a56767174
# ╠═4d8089c3-435e-4b1a-a185-bd79e1425c13
# ╠═9fc7b12c-b770-4431-bc38-0293b88a13f1
# ╠═7e4c3269-2a2b-4b15-a373-4daf815a72ef
# ╟─8b13e40c-b2fd-4f2d-856e-7605c5f5c0b7
# ╟─f00fb7d4-5183-483a-90a2-73ff9be3a24b
# ╟─efd98ea7-b86a-41e1-9eca-6675dffc7edd
# ╠═c3837d06-7c00-49ca-8592-6b9e8e7ef528
# ╠═1870556b-2e19-4e64-9df1-338daae29805
# ╠═2d3b75aa-da53-4490-a743-9a6ca405759f
# ╠═ade72d79-4082-4c4e-b14a-744651731d9d
# ╟─320405f5-4feb-4c61-98f4-ff13643490ae
# ╟─122317ec-7c00-4b1e-9045-1314f444e56d
# ╟─017fe43a-df10-40d9-b9d1-d93d71500bd0
# ╠═d37b1618-f0ed-42af-a336-8f9a341f2445
# ╠═7635fbef-61fa-4311-bda4-f7f384a4832d
# ╟─11c49d54-2088-4f6c-83a1-129960889492
# ╟─a91fd5a5-cd7b-4676-a6c3-13d8c81453d9
# ╟─c624e494-ade4-4a05-98e7-f7beaeaf847f
# ╠═a0609e90-45ea-4494-ad0a-86d339dff917
# ╟─581e9601-c169-4b2a-acec-391a75be26da
# ╠═15a6ff92-a3c1-4ff1-bc4b-674adbc53c0d
# ╟─1e763d9c-44e8-4614-94fb-be9b65f100a9
# ╟─7484f5a6-0cbc-451f-b790-aea05813c425
# ╟─fac54db6-e2a1-4956-8b7c-091938e74b8e
# ╠═29b04989-caa7-4e73-a9a1-cdad9e3509d8
# ╟─982837e3-6d4c-4012-a014-886d38d19121
# ╟─107e9699-6baa-4978-a68f-9fc22d76d6a1
# ╟─ac9da8e2-aee9-4c79-8e7d-f35eb7e674f2
# ╟─4a14052c-5327-43a4-bf2e-ad739bf5213d
# ╟─140d7064-df9e-4dc8-886e-edb381792164
# ╟─78d99536-b253-491e-b08d-df7cb97a2063
# ╟─7778b027-68be-46c7-a8f3-31413fa5554f

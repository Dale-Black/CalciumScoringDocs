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

# ╔═╡ 359ca243-8ad1-4728-97fe-35022314dae3
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
		Pkg.add("Statistics")
		Pkg.add("DICOM")
		Pkg.add("Unitful")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI
	using CairoMakie
	using CalciumScoring
	
	using CSV: read
	using DataFrames: DataFrame
	using ImageMorphology: dilate
	using Statistics: mean
	using DICOM: dcmdir_parse, @tag_str
	using Unitful: mm, mg, ustrip

	TableOfContents()
end

# ╔═╡ f5fd4bdc-1093-11ee-35c9-db1ff5709044
md"""
# Overview
[Previously](01-getting-started.jl), we introduced the CalciumScoring.jl package and specifically showcased the Agatston scoring method. In this notebook we will examine how to use the Volume Fraction Calium Mass method. To understand the theory behind this technique, please see [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1101/2023.01.12.23284482)
"""

# ╔═╡ d53b940a-978a-4497-a274-85cc31d6d12a
md"""
## Preparation
Let's quickly (1) import all the neccessary packages, (2) load the same simulated phantom images and masks from before, and (3) calculate the ground truth calcium mass (using Unitful to keep track of units).
"""

# ╔═╡ 3ce4d714-6852-4f73-b5c2-0d67b6405d9a
md"""
# Volume Fraction Calcium Mass
This calcium quantification technique does not require any intensity-based threhsolding. To calculate the calclium contained within a region of interest (ROI), the mean intensity of the background material is needed.

To find the background intensity, we will first create a background ring mask by dilating the dilated mask again and then subtracting.
"""

# ╔═╡ 61a0988b-6776-4620-83a5-4f274aa73952
md"""
## Background calculation
"""

# ╔═╡ 9d266d02-cd20-4db1-83b6-d18ca5a68006
md"""
## Visualize Insert and Background Masks
"""

# ╔═╡ 7f612715-b4fc-4de3-be11-ddc6c2ed72c1
md"""
## Results
"""

# ╔═╡ fec583af-b4a0-452d-9ac5-5b953caef70c
md"""
# Next Steps
We just demonstrated how `score()` can be used with the `VolumeFraction()` algorithm (along with units from Unitful.jl). This is the most well tested calcium scoring algorithm in the library, but check out the [Integrated Calcium Mass](\03-integated.jl) tutorial to see how to implement another approach.
"""

# ╔═╡ c916adfb-9336-42cc-ae6a-9b1f92e894a2
md"""
#### Helper Functions
"""

# ╔═╡ 8b9a3224-5306-4648-b51e-5792ce9e0e55
function dilate_recursively(mask, n)
    dilated_mask = copy(mask)
    for _ in 1:n
        dilated_mask = dilate(dilated_mask)
    end
    return dilated_mask
end

# ╔═╡ a084d53a-cb7f-4d65-a716-7d35a4b47b6f
begin
	# Validation Phantom
	dcm_path = joinpath(pwd(), "dcms/val/52_59_73/medium/80");
	dcms = dcmdir_parse(dcm_path)
	dcm_array = cat([dcms[i][(0x7fe0, 0x0010)] for i in 1:length(dcms)]...; dims=3)

	# Calcification Masks
	mask_large_insert_high_density = Array(
		read(
			joinpath(pwd(), "dcms", "val_masks", "medium", "mask_L_LD.csv"), DataFrame;
			header=false
		)
	)
	mask_large_insert_high_density_3D = cat(
		mask_large_insert_high_density,
		mask_large_insert_high_density,
		mask_large_insert_high_density; 
		dims = 3
	)

	# Calibration Rod
	dcm_path_calibration = joinpath(pwd(), "dcms/cal/200/30/80");
	dcms_calibration = dcmdir_parse(dcm_path_calibration)
	dcm_array_calibration = cat(
		[dcms_calibration[i][(0x7fe0, 0x0010)] for i in 1:length(dcms_calibration)]...;
		dims=3
	)

	center_insert1, center_insert2 = 187, 318
	calibration_rod = zeros(25, 25, size(dcm_array_calibration, 3))
				
	for z in axes(dcm_array_calibration, 3)
		rows, cols, depth = size(dcm_array_calibration)
		half_row, half_col = center_insert1, center_insert2
		offset = 12
		row_range = half_row-offset:half_row+offset
		col_range = half_col-offset:half_col+offset	
		calibration_rod[:, :, z] .= dcm_array_calibration[row_range, col_range, z];
	end
	hu_calcium, ρ_calcium = mean(calibration_rod), 0.200mg/mm^3


	# Dilate mask
	dilated_mask_large_high_density = dilate_recursively(mask_large_insert_high_density, 2)
	dilated_mask_large_high_density_3D = cat(
		dilated_mask_large_high_density, 
		dilated_mask_large_high_density,
		dilated_mask_large_high_density;
		dims=3)
end;

# ╔═╡ 70f4abc8-36a3-4f2c-9e54-21562c87d16f
begin
	slice_thickness = dcms[1].meta[tag"Slice Thickness"]mm
	pixel_spacing = dcms[1][tag"Pixel Spacing"]mm
	voxel_spacing = [pixel_spacing..., slice_thickness]

	num_slices = length(axes(dcm_array, 3))
	ρ_high_density = 0.073mg/mm^3
	radius = (5/2)mm
	calcium_volume_large = π * radius^2 * slice_thickness * num_slices
	calcium_mass_large = calcium_volume_large * ρ_high_density
end

# ╔═╡ 3abad7b7-37b4-415e-b3fe-34428d2b8d17
voxel_size = voxel_spacing[1] * voxel_spacing[2] * voxel_spacing[3]

# ╔═╡ a0a67d7a-183f-464e-b21a-032d32de09bb
@bind z1 PlutoUI.Slider(axes(dcm_array, 3); default = 2, show_value = true)

# ╔═╡ e43b6288-08f8-463a-8324-89eb32bcdaba
begin
	_background_mask_3D = dilate_recursively(dilated_mask_large_high_density_3D, 3)
	background_mask_3D = Bool.(_background_mask_3D - dilated_mask_large_high_density_3D)
end;

# ╔═╡ c37667c7-5074-43a1-a75d-8e44884c7160
let
	idxs = getindex.(findall(isone, dilated_mask_large_high_density_3D[:, :, z1]), [1 2])
	idxs2 = getindex.(findall(isone, background_mask_3D[:, :, z1]), [1 2])
	
	f = Figure(resolution = (1200, 800))
	ax = Axis(
		f[1, 1],
		title = "DICOM Array"
	)
	heatmap!(transpose(dcm_array[:, :, z1]); colormap = :grays)
	scatter!(idxs[:, 2], idxs[:, 1]; markersize = 4, color = (:red, 0.2))
	hidedecorations!(ax)
	scatter!(idxs2[:, 2], idxs2[:, 1]; markersize = 4, color = (:blue, 0.2))
	hidedecorations!(ax)
	f
end

# ╔═╡ 539dfd12-5b93-44d4-b883-0d08ebadbe55
hu_background = mean(dcm_array[background_mask_3D])

# ╔═╡ 5be9d073-e8c5-44d6-8982-639cbe8c3c6a
volume_fraction_mass = score(dcm_array[dilated_mask_large_high_density_3D], hu_calcium, hu_background, voxel_size, ρ_calcium, VolumeFraction())

# ╔═╡ 050f4296-e9ea-4bc7-ae28-b98f2a6c0211
md"""
Let's compare that with the ground truth mass. We see that the ground truth mass = $(round.(mg, calcium_mass_large; digits = 2)) is close to the calculated mass = $(round.(mg, volume_fraction_mass; digits = 2)). This is a better estimate than the previous [Agatston scoring](01-getting-started.jl) calculation.
"""

# ╔═╡ Cell order:
# ╟─f5fd4bdc-1093-11ee-35c9-db1ff5709044
# ╟─d53b940a-978a-4497-a274-85cc31d6d12a
# ╠═359ca243-8ad1-4728-97fe-35022314dae3
# ╠═a084d53a-cb7f-4d65-a716-7d35a4b47b6f
# ╠═70f4abc8-36a3-4f2c-9e54-21562c87d16f
# ╟─3ce4d714-6852-4f73-b5c2-0d67b6405d9a
# ╟─61a0988b-6776-4620-83a5-4f274aa73952
# ╠═e43b6288-08f8-463a-8324-89eb32bcdaba
# ╟─9d266d02-cd20-4db1-83b6-d18ca5a68006
# ╟─a0a67d7a-183f-464e-b21a-032d32de09bb
# ╟─c37667c7-5074-43a1-a75d-8e44884c7160
# ╠═539dfd12-5b93-44d4-b883-0d08ebadbe55
# ╟─7f612715-b4fc-4de3-be11-ddc6c2ed72c1
# ╠═3abad7b7-37b4-415e-b3fe-34428d2b8d17
# ╠═5be9d073-e8c5-44d6-8982-639cbe8c3c6a
# ╟─050f4296-e9ea-4bc7-ae28-b98f2a6c0211
# ╟─fec583af-b4a0-452d-9ac5-5b953caef70c
# ╟─c916adfb-9336-42cc-ae6a-9b1f92e894a2
# ╟─8b9a3224-5306-4648-b51e-5792ce9e0e55

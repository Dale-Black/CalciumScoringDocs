### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 9374f51f-6b55-4c91-b063-1b96e4d17fbd
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

# ╔═╡ 360d6b4b-da76-43b3-a6a0-96f51e14b224
md"""
# Overview
[Previously](02-volume-fraction.jl), we looked into the Volume Fraction Calcium Mass quantification approach. In this notebook we will examine how to use the Integrated Calium Mass method. To understand the theory behind this technique, please see [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1101/2023.01.12.23284482)
and [Integrated intensity-based technique for coronary artery calcium mass measurement: A phantom study](https://doi.org/10.1002/mp.16326)
"""

# ╔═╡ f7b41159-befc-4806-9650-97a122e31428
md"""
## Preparation
Let's quickly (1) import all the neccessary packages, (2) load the same simulated phantom images and masks from before, and (3) calculate the ground truth calcium mass (using Unitful to keep track of units).
"""

# ╔═╡ e0b4c970-09da-4b19-baab-6957e7932336
md"""
#### Helper Functions
"""

# ╔═╡ 1f30753c-6a62-4e78-b9c0-3b6a5fd4e3bb
function dilate_recursively(mask, n)
    dilated_mask = copy(mask)
    for _ in 1:n
        dilated_mask = dilate(dilated_mask)
    end
    return dilated_mask
end

# ╔═╡ 82f5c527-c5fd-4f75-8f43-16afa68ab9af
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

# ╔═╡ f4b2034e-913f-4dbd-86e2-c1bc5f2b3d79
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

# ╔═╡ Cell order:
# ╟─360d6b4b-da76-43b3-a6a0-96f51e14b224
# ╟─f7b41159-befc-4806-9650-97a122e31428
# ╠═9374f51f-6b55-4c91-b063-1b96e4d17fbd
# ╠═82f5c527-c5fd-4f75-8f43-16afa68ab9af
# ╠═f4b2034e-913f-4dbd-86e2-c1bc5f2b3d79
# ╟─e0b4c970-09da-4b19-baab-6957e7932336
# ╟─1f30753c-6a62-4e78-b9c0-3b6a5fd4e3bb

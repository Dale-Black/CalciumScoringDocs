### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ f3363153-9287-482e-8049-bbc317aa8b51
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
		Pkg.add("ImageFiltering")
		Pkg.add("Noise")
		Pkg.add("Unitful")
		Pkg.add(url="https://github.com/Dale-Black/CalciumScoring.jl")
	end
	
	using PlutoUI, CSV, DataFrames, CairoMakie, ImageMorphology, GLM, Statistics, ImageFiltering, Noise, Unitful, CalciumScoring
	using Unitful: mm, mg
end

# ╔═╡ 0a38e736-a860-49b6-9152-a51f33012117
md"""
# Overview
This tutorial will show how to calculate the calcium mass in a non-contrast dual-energy CT image using the `MaterialDecomposition()` technique from CalciumScoring.jl
"""

# ╔═╡ 805556cc-0ed5-4b60-909b-933626c93e95
md"""
## Import Packages
First, let's import the most up-to-date version of CalciumScoring.jl, which can be found on the main/master branch of the [GitHub repository](https://github.com/Dale-Black/CalciumScoring.jl).

Because we are using the unregistered version (most recent), we will need to `Pkg.add` this explicitly without using Pluto's built-in package manager. Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.

To help with the formatting of this documentation, we will also add [PlutoUI.jl](https://github.com/JuliaPluto/PlutoUI.jl). We will also create some simulated coronary artery images, which will require [ImageFiltering.jl](https://github.com/JuliaImages/ImageFiltering.jl) and [Noise.jl](https://github.com/roflmaostc/Noise.jl) Lastly, to visualize the results, we're going to add [Makie.jl](https://github.com/JuliaPlots/Makie.jl).
"""

# ╔═╡ a77307a4-7d8e-4e40-8757-952a34d4d423
TableOfContents()

# ╔═╡ d74d15b9-b330-4f31-af47-88a7f8d314b1
md"""
# Calculate Ground Truth Values

Below, in the appendix, we created simulated dual-energy images. One set of calibration images (high and low energy) and another set of measurement images (high and low energy).

Since we created the calcifications as a circular mask, we know the exact geometry of the ideal calcium mass. We can multiply this by the voxel size to determine the ground truth calcium volume. Lastly, we can take the density of calcium and multiply it by the volume to determine the mass of the calcium
"""

# ╔═╡ 2f33b7b1-1a52-41e3-8521-f15f5594e971
md"""
# Determine Calibration Parameters
To determine the calibration parameters for calcium mass measurement, we first require a simulated calibration dual-energy image (low and high energy), with various calcium rods and known calcium densities. Below, in the appendix, we created such images. Let's quickly visualize these images below.
"""

# ╔═╡ a528497a-0b2d-4561-836f-69a57d94e0db
md"""
## Measure Attenuations
We will use the masks to measure the calibration rod intensities for both low and high-energy images. We want to avoid any partial volume effect, so we will erode the masks before measuring the mean intensities.
"""

# ╔═╡ 8dbf0852-8139-45d1-acff-246202396715
md"""
## Find Parameters
These calculated mean intensities can then be used to fit the parameters in equation below

```math
\begin{aligned}	F = \frac{a_o + a_1x + a_2y + a_3x^2 + a_4xy + a_5y^2}{1 + b_1x + b_2y} 
\end{aligned}
\tag{1}
```
"""

# ╔═╡ 239a4e84-48a8-48b2-a45b-32ed9d512ad5
calcium_densities = [0, 0.025, 0.1, 0.25, 0.35, 0.45, 0.60, 0.80]

# ╔═╡ 25faf325-2ccb-47b2-80f8-a44882d53d8a
md"""
## Check Fit
To make sure the equation was properly fit, we can double-check that the calculated densities on the calibration rods are nearly identical to the ground truth densities when utilizing the parameters above.
"""

# ╔═╡ b6cdb1e7-be55-4898-91d2-624e6901965e
md"""
# Measure Calcium
Now that the given parameters seem to work, we can measure an unknown calcification using the given parameters. First, we want to mask the calcification ROI and ensure the ROI includes all voxels containing any calcium. For this, we will dilate the calcification mask.
"""

# ╔═╡ 5a720fc1-d082-4f0b-b7ea-fbcc6a978345
md"""
## Measure Attenuations
Now, we will measure the average attenuation of the dilated mask on the low and high-energy calcifications. These intensity measurements and the above parameters will allow us to calculate the density of the given calcification ROI. Assuming we know the size of the voxels, we can then convert this relative density to an absolute mass by multiplying the calculated density by the volume of the ROI.
"""

# ╔═╡ 4384ae7f-ed4a-47ec-9979-1a7a25409fe7
md"""
## Results
"""

# ╔═╡ 5d92ebf6-ef70-42d6-ab4b-3ce5d516f4fe
md"""
# Appendix
"""

# ╔═╡ 296b5f5f-e8db-45a9-8010-3b98056d6781
md"""
## Create Simulated Calibration Images
To understand calcium scoring, we will create simulated CT coronary artery calcium images.

First, let's create a 2D coronary artery calcification with a density of 0.2 ``mg/mm^{3}`` and an intensity value of 300 HU. We will place this in simulated heart tissue with a simulated intensity value of 40 HU. We will also prepare each voxel to be of size = [0.5, 0.5, 0.5] ``mm^{3}``
"""

# ╔═╡ 11e5c831-a033-4a7c-8019-8daebf5e92c3
begin
	# ρ_heart_tissue = 1.055mg/mm^3
	ρ_calcium = 0.2mg/mm^3

	hu_heart_tissue_low = 40 # HU
	hu_heart_tissue_high = 32 # HU
	hu_calcium_low = 376 # HU
	hu_calcium_high = 282 # HU

	spacing = [0.5, 0.5, 0.5]mm
	voxel_size = spacing[1] * spacing[2] * spacing[3]
end

# ╔═╡ 25c56212-e9e4-46fc-8cdd-f8d80fc7661d
function create_circular_mask(h, w, center_circle, radius_circle)
    Y, X = collect(1:h), collect(1:w)'
    dist_from_center = sqrt.((X .- center_circle[1]) .^ 2 .+ (Y .- center_circle[2]) .^ 2)
    mask = dist_from_center .<= radius_circle
    return mask
end

# ╔═╡ 99e4b65f-c8bf-4a60-89d0-2f3d036a05dd
begin
	h, w, rad = 400, 400, 15

	mask0 = create_circular_mask(h, w, (300, 50), rad)

	mask1 = create_circular_mask(h, w, (50, 50), rad)
	rod1_low = mask1 * 69 # 25 mg/cc
	rod1_high = mask1 * 55 # 25 mg/cc
	
	mask2 = create_circular_mask(h, w, (100, 100), rad)
	rod2_low = mask2 * 201 # 100 mg/cc
	rod2_high = mask2 * 152 # 100 mg/cc
	
	mask3 = create_circular_mask(h, w, (150, 150), rad)
	rod3_low = mask3 * 463 # 250 mg/cc
	rod3_high = mask3 * 345 # 250 mg/cc
	
	mask4 = create_circular_mask(h, w, (200, 200), rad)
	rod4_low = mask4 * 633 # 350 mg/cc
	rod4_high = mask4 * 471 # 350 mg/cc
	
	mask5 = create_circular_mask(h, w, (250, 250), rad)
	rod5_low = mask5 * 803 # 450 mg/cc
	rod5_high = mask5 * 596 # 450 mg/cc
	
	mask6 = create_circular_mask(h, w, (300, 300), rad)
	rod6_low = mask6 * 1050 # 600 mg/cc
	rod6_high = mask6 * 780 # 600 mg/cc
		
	mask7 = create_circular_mask(h, w, (350, 350), rad)
	rod7_low = mask7 * 1379 # 800 mg/cc
	rod7_high = mask7 * 1022 # 800 mg/cc

	rods_low = rod1_low + rod2_low + rod3_low + rod4_low + rod5_low + rod6_low + rod7_low
	rods_high = rod1_high + rod2_high + rod3_high + rod4_high + rod5_high + rod6_high + rod7_high
end;

# ╔═╡ 744b47ef-2d5f-4f7b-989f-465d9ac7fc40
begin	
	mask0_3D = cat(mask0, mask0, dims=3)
	mask0_3D = cat(mask0_3D, mask0_3D, dims=3)
	
	mask1_3D = cat(mask1, mask1, dims=3)
	mask1_3D = cat(mask1_3D, mask1_3D, dims=3)

	mask2_3D = cat(mask2, mask2, dims=3)
	mask2_3D = cat(mask2_3D, mask2_3D, dims=3)

	mask3_3D = cat(mask3, mask3, dims=3)
	mask3_3D = cat(mask3_3D, mask3_3D, dims=3)

	mask4_3D = cat(mask4, mask4, dims=3)
	mask4_3D = cat(mask4_3D, mask4_3D, dims=3)

	mask5_3D = cat(mask5, mask5, dims=3)
	mask5_3D = cat(mask5_3D, mask5_3D, dims=3)

	mask6_3D = cat(mask6, mask6, dims=3)
	mask6_3D = cat(mask6_3D, mask6_3D, dims=3)

	mask7_3D = cat(mask7, mask7, dims=3)
	mask7_3D = cat(mask7_3D, mask7_3D, dims=3)
end;

# ╔═╡ 0fb61ccb-43b7-4b1e-8b92-9692d0a42548
begin
	pure_calcification_low = @. ifelse(rods_low == 0, rods_low .+ hu_heart_tissue_low, rods_low)
	pure_calcification_high = @. ifelse(rods_low == 0, rods_high .+ hu_heart_tissue_high, rods_high)
end;

# ╔═╡ eb39e77a-1256-4006-a625-b9dc4bd114eb
md"""
Now, let's simulate 3D by stacking this slice along the third dimension
"""

# ╔═╡ 6c22fc99-923a-486f-9d34-8c6d04f5e969
begin
	pure_calcification_3D_low = cat(pure_calcification_low, pure_calcification_low, pure_calcification_low, pure_calcification_low, dims=3)
	pure_calcification_3D_high = cat(pure_calcification_high, pure_calcification_high, pure_calcification_high, pure_calcification_high, dims=3)
end;

# ╔═╡ 6e652ad4-de3c-47cd-b8c6-9b0c2ebacee3
md"""
Finally, let's add a gaussian filter and poisson noise to simulate a more realistic coronary artery calcification
"""

# ╔═╡ 1dd08bcb-e945-4e9c-9bbe-3fc37a843894
begin
	gaussian_calcification1_low = imfilter(pure_calcification_low, Kernel.gaussian(3))
	gaussian_calcification2_low = imfilter(pure_calcification_low, Kernel.gaussian(3))
	gaussian_calcification3_low = imfilter(pure_calcification_low, Kernel.gaussian(3))
	gaussian_calcification4_low = imfilter(pure_calcification_low, Kernel.gaussian(3))
	
	gaussian_poisson_calcification1_low = poisson(gaussian_calcification1_low, 100)
	gaussian_poisson_calcification2_low = poisson(gaussian_calcification2_low, 100)
	gaussian_poisson_calcification3_low = poisson(gaussian_calcification3_low, 100)
	gaussian_poisson_calcification4_low = poisson(gaussian_calcification4_low, 100)

	
	gaussian_poisson_calcification_3D_low = cat(gaussian_poisson_calcification1_low, gaussian_poisson_calcification2_low, gaussian_poisson_calcification3_low, gaussian_poisson_calcification4_low, dims=3)
end;

# ╔═╡ b8977e39-f5c4-4d3d-a6e2-104d5f11b549
@bind x PlutoUI.Slider(axes(gaussian_poisson_calcification_3D_low, 3); default=2, show_value=true)

# ╔═╡ 48ca6e64-010f-4460-8a2f-2eaf06cde1ec
means1 = [
	mean(gaussian_poisson_calcification_3D_low[erode(erode(mask0_3D))]), mean(gaussian_poisson_calcification_3D_low[erode(erode(mask1_3D))]), mean(gaussian_poisson_calcification_3D_low[erode(erode(mask2_3D))]), mean(gaussian_poisson_calcification_3D_low[erode(erode(mask3_3D))]), mean(gaussian_poisson_calcification_3D_low[erode(erode(mask4_3D))]), mean(gaussian_poisson_calcification_3D_low[erode(erode(mask5_3D))]), mean(gaussian_poisson_calcification_3D_low[erode(erode(mask6_3D))]),
	mean(gaussian_poisson_calcification_3D_low[erode(erode(mask7_3D))])
]

# ╔═╡ 3ab8f65a-ae68-49b0-b717-32de09d4b56f
begin
	gaussian_calcification1_high = imfilter(pure_calcification_high, Kernel.gaussian(3))
	gaussian_calcification2_high = imfilter(pure_calcification_high, Kernel.gaussian(3))
	gaussian_calcification3_high = imfilter(pure_calcification_high, Kernel.gaussian(3))
	gaussian_calcification4_high = imfilter(pure_calcification_high, Kernel.gaussian(3))
	
	gaussian_poisson_calcification1_high = poisson(gaussian_calcification1_high, 100)
	gaussian_poisson_calcification2_high = poisson(gaussian_calcification2_high, 100)
	gaussian_poisson_calcification3_high = poisson(gaussian_calcification3_high, 100)
	gaussian_poisson_calcification4_high = poisson(gaussian_calcification4_high, 100)

	
	gaussian_poisson_calcification_3D_high = cat(gaussian_poisson_calcification1_high, gaussian_poisson_calcification2_high, gaussian_poisson_calcification3_high, gaussian_poisson_calcification4_high, dims=3)
end;

# ╔═╡ e27cf9d6-62b1-49d6-9bc1-07f7f6711af0
means2 = [
	mean(gaussian_poisson_calcification_3D_high[erode(erode(mask0_3D))]), mean(gaussian_poisson_calcification_3D_high[erode(erode(mask1_3D))]), mean(gaussian_poisson_calcification_3D_high[erode(erode(mask2_3D))]), mean(gaussian_poisson_calcification_3D_high[erode(erode(mask3_3D))]), mean(gaussian_poisson_calcification_3D_high[erode(erode(mask4_3D))]), mean(gaussian_poisson_calcification_3D_high[erode(erode(mask5_3D))]), mean(gaussian_poisson_calcification_3D_high[erode(erode(mask6_3D))]),
	mean(gaussian_poisson_calcification_3D_high[erode(erode(mask7_3D))])
]

# ╔═╡ f5ca6a29-bbbb-4659-ab2e-ea93a5dfdafc
calculated_intensities = hcat(means1, means2) # low energy, high energy

# ╔═╡ 9cfbd6b5-3db1-489d-b57f-e7e7a48a0068
ps = fit_calibration(calculated_intensities, calcium_densities)

# ╔═╡ dc58098b-8756-4119-9974-9335b5aeb5da
begin
	predicted_densities = []
	
	for i in 1:length(calcium_densities)
		append!(
			predicted_densities, 
			score(calculated_intensities[i, 1], calculated_intensities[i, 2], ps, MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ cd417d4c-2c08-47d8-9925-f0eeaf063f63
df = DataFrame(
	calcium_densities = calcium_densities,
	predicted_densities = predicted_densities,
	mean_intensities_low = means1,
	mean_intensities_high = means2,
)

# ╔═╡ 4c0bc41a-20bf-497d-8dbc-3b12389cb1aa
md"""
### Visualize
Below we will visualize the pure coronary artery calcification along with a more realistic simulated coronary artery calcification with added noise and blurring
"""

# ╔═╡ 321f8a2a-224e-48e2-9fb8-c8a8922c3e4a
md"""
#### Coronary Artery Calcification (Ideal)
"""

# ╔═╡ 02721b3d-3a23-4633-9706-40c94856b91d
@bind a PlutoUI.Slider(axes(pure_calcification_3D_low, 3); default=2, show_value=true)

# ╔═╡ 9cfd3199-e014-4c38-97d4-a2493e0ac437
let
	f = Figure(resolution=(800, 500))
	ax = CairoMakie.Axis(
		f[1, 1],
		title = "Low Energy Calibration Image"
	)
	heatmap!(gaussian_poisson_calcification_3D_low[:, :, a], colormap=:grays)

	ax = CairoMakie.Axis(
		f[1, 2],
		title = "High Energy Calibration Image"
	)
	heatmap!(gaussian_poisson_calcification_3D_high[:, :, a], colormap=:grays)

	f
end

# ╔═╡ c7e6d66d-18f9-4795-89ed-11c5c312d8f4
md"""
#### Low energy
"""

# ╔═╡ c97c4a72-7857-42f7-abe4-8ad247ff538c
heatmap(pure_calcification_3D_low[:, :, a], colormap=:grays)

# ╔═╡ b3d6ec94-4cc2-43d4-b64f-f4d6bf12d54d
md"""
#### High energy
"""

# ╔═╡ ee56071a-0d3d-48c3-9aff-73933c0c0146
heatmap(pure_calcification_3D_high[:, :, a], colormap=:grays)

# ╔═╡ 167385d2-8d62-4193-885e-4bbc3683461d
md"""
#### Coronary Artery Calcification (Noisy)
"""

# ╔═╡ 6d6076eb-e1bc-44f3-a913-cda403af0d70
@bind b1 PlutoUI.Slider(axes(pure_calcification_3D_low, 3); default=2, show_value=true)

# ╔═╡ 4f031d7f-70c1-42d9-af64-8057178b7fca
md"""
#### Low energy
"""

# ╔═╡ 52acce68-1f6e-426b-ab84-c52e60b6870b
heatmap(gaussian_poisson_calcification_3D_low[:, :, b1], colormap=:grays)

# ╔═╡ 66d311a5-0dd9-4d4e-99fd-377c3bfab82f
md"""
#### High energy
"""

# ╔═╡ cc9000bb-dc0c-4cdd-b15a-f2af91da5ea8
heatmap(gaussian_poisson_calcification_3D_high[:, :, b1], colormap=:grays)

# ╔═╡ f9c17967-6f97-44bb-9c8f-dbc35f09b0e2
md"""
## Create Simulated Measurement Image
Once the above calibration is computed for a given CT system, the parameters are all that is needed for calcium mass measurement.
"""

# ╔═╡ 0a7c1547-beba-4562-9cd7-c2f882612182
mask = create_circular_mask(h, w, (200, 200), rad);

# ╔═╡ f394d548-3792-4999-bb54-5208e46e6a34
begin
	mask_3D = cat(mask, mask, dims=3)
	mask_3D = cat(mask_3D, mask_3D, dims=3)
end;

# ╔═╡ 080b868a-b8c2-445c-b53f-d587ac6e3288
ground_truth_number_voxels = length(findall(x -> x == 1, mask_3D))

# ╔═╡ 04b354df-36cd-422c-b5a1-09f19a0662ed
ground_truth_calcium_volume = ground_truth_number_voxels * voxel_size

# ╔═╡ dd93a5d6-fa5f-4a80-9b76-79abe755e833
ground_truth_calcium_mass = ground_truth_calcium_volume * ρ_calcium

# ╔═╡ a0cfcbd7-0f49-472d-a6e2-2de7b87bfb5e
dilated_mask_3D = dilate(dilate(dilate(dilate(dilate(dilate(dilate(mask_3D)))))));

# ╔═╡ 8b0dd679-6266-4ce7-960b-db56ddef1afa
vol_ROI = sum(dilated_mask_3D) * spacing[1] * spacing[2] * spacing[3]

# ╔═╡ b3a78f50-d709-40b8-a5f4-840534939091
begin
	pure_calcification_low_measure = @. ifelse(mask == 0, mask + hu_heart_tissue_low, mask + hu_calcium_low)
	pure_calcification_high_measure = @. ifelse(mask == 0, mask + hu_heart_tissue_high, mask + hu_calcium_high)
end;

# ╔═╡ fa77e270-ab83-4e2f-a304-9f2ca7e44c0c
md"""
Now, let's simulate 3D by stacking this slice along the third dimension
"""

# ╔═╡ a323a639-d40f-4f68-8e2b-822e6d70177c
begin
	pure_calcification_3D_low_measure = cat(pure_calcification_low_measure, pure_calcification_low_measure, pure_calcification_low_measure, pure_calcification_low_measure, dims=3)
	
	pure_calcification_3D_high_measure = cat(pure_calcification_high_measure, pure_calcification_high_measure, pure_calcification_high_measure, pure_calcification_high_measure, dims=3)
end;

# ╔═╡ aaa53ab7-2aca-46c5-b15a-4506ec061ebd
md"""
Finally, let's add a gaussian filter and poisson noise to simulate a more realistic coronary artery calcification
"""

# ╔═╡ 0a8f0a44-33fd-470f-b5fb-ca0138b21b5c
begin
	gaussian_calcification1_low_measure = imfilter(pure_calcification_low_measure, Kernel.gaussian(3))
	gaussian_calcification2_low_measure = imfilter(pure_calcification_low_measure, Kernel.gaussian(3))
	gaussian_calcification3_low_measure = imfilter(pure_calcification_low_measure, Kernel.gaussian(3))
	gaussian_calcification4_low_measure = imfilter(pure_calcification_low_measure, Kernel.gaussian(3))
	
	gaussian_poisson_calcification1_low_measure = poisson(gaussian_calcification1_low_measure, 100)
	gaussian_poisson_calcification2_low_measure = poisson(gaussian_calcification2_low_measure, 100)
	gaussian_poisson_calcification3_low_measure = poisson(gaussian_calcification3_low_measure, 100)
	gaussian_poisson_calcification4_low_measure = poisson(gaussian_calcification4_low_measure, 100)

	
	gaussian_poisson_calcification_3D_low_measure = cat(gaussian_poisson_calcification1_low_measure, gaussian_poisson_calcification2_low_measure, gaussian_poisson_calcification3_low_measure, gaussian_poisson_calcification4_low_measure, dims=3)
end;

# ╔═╡ 964b6d60-2e1b-468f-9d5c-7ce2c631e8c7
low_energy_intensity = mean(gaussian_poisson_calcification_3D_low_measure[dilated_mask_3D])

# ╔═╡ d039d82e-3f26-4434-a9f8-cabad86b8771
begin
	gaussian_calcification1_high_measure = imfilter(pure_calcification_high_measure, Kernel.gaussian(3))
	gaussian_calcification2_high_measure = imfilter(pure_calcification_high_measure, Kernel.gaussian(3))
	gaussian_calcification3_high_measure = imfilter(pure_calcification_high_measure, Kernel.gaussian(3))
	gaussian_calcification4_high_measure = imfilter(pure_calcification_high_measure, Kernel.gaussian(3))
	
	gaussian_poisson_calcification1_high_measure = poisson(gaussian_calcification1_high_measure, 100)
	gaussian_poisson_calcification2_high_measure = poisson(gaussian_calcification2_high_measure, 100)
	gaussian_poisson_calcification3_high_measure = poisson(gaussian_calcification3_high_measure, 100)
	gaussian_poisson_calcification4_high_measure = poisson(gaussian_calcification4_high_measure, 100)

	
	gaussian_poisson_calcification_3D_high_measure = cat(gaussian_poisson_calcification1_high_measure, gaussian_poisson_calcification2_high_measure, gaussian_poisson_calcification3_high_measure, gaussian_poisson_calcification4_high_measure, dims=3)
end;

# ╔═╡ 68f05736-722c-49b8-9222-8c5548599cbf
high_energy_intensity = mean(gaussian_poisson_calcification_3D_high_measure[dilated_mask_3D])

# ╔═╡ b4269f46-bdb4-4da6-acce-2f31f143d36a
begin
	material_decomposition_mass = score(low_energy_intensity, high_energy_intensity, ps, vol_ROI, MaterialDecomposition())
	material_decomposition_mass =  (material_decomposition_mass)mg/mm^3
end

# ╔═╡ 87c3f233-3607-4805-9394-16b6857a5887
md"""
Let's compare that with the ground truth mass

We see that the ground truth mass = $(ground_truth_calcium_mass) is close to the calculated mass = $(material_decomposition_mass)
"""

# ╔═╡ a81771c1-78ed-40b6-8517-ec71c5cb17a7
md"""
### Visualize
Below we will visualize the pure coronary artery calcification along with a more realistic simulated coronary artery calcification with added noise and blurring
"""

# ╔═╡ 6e33a03b-6484-49de-833c-7337b2cb1e2d
md"""
#### Coronary Artery Calcification (Ideal)
"""

# ╔═╡ 939ef2cc-11b2-40ae-b489-59d15fe5a437
@bind a_measure PlutoUI.Slider(axes(pure_calcification_3D_low_measure, 3); default=2, show_value=true)

# ╔═╡ de246ed9-912d-4c55-bd71-09972abb3d02
md"""
#### Low Energy
"""

# ╔═╡ 653b3f92-e48e-4440-8913-fca7e90ca3fc
heatmap(pure_calcification_3D_low_measure[:, :, a_measure], colormap=:grays)

# ╔═╡ 2106147b-e6c6-48fd-8e65-ce1b329c06ae
md"""
#### High Energy
"""

# ╔═╡ 91b745c3-689f-4164-b3e5-4f49eecb93a4
heatmap(pure_calcification_3D_high_measure[:, :, a_measure], colormap=:grays)

# ╔═╡ 4a64c4c4-f713-4eb1-a18d-809e9816d594
md"""
#### Coronary Artery Calcification (Noisy)
"""

# ╔═╡ f026dcbf-0e74-49a0-94c0-32cfd77b869c
@bind b_measure PlutoUI.Slider(axes(gaussian_poisson_calcification_3D_low_measure, 3); default=2, show_value=true)

# ╔═╡ bc224643-f914-42ee-96ea-ce5f40472350
md"""
#### Low Energy
"""

# ╔═╡ 1add336d-dbcf-4466-a25a-043bde36a021
heatmap(gaussian_poisson_calcification_3D_low_measure[:, :, b_measure], colormap=:grays)

# ╔═╡ 4d8d6448-b60c-4072-959f-d32310c5a9dd
md"""
#### High Energy
"""

# ╔═╡ a1b37552-936a-4d56-ba16-d280b29adaa8
heatmap(gaussian_poisson_calcification_3D_high_measure[:, :, b_measure], colormap=:grays)

# ╔═╡ f95221cd-ae8d-4cf0-bfe9-cb7e005b242c
md"""
#### Helper functions
"""

# ╔═╡ 40ef54ea-3afd-487b-8bea-806e33134b68
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 147b5c20-b93d-45d8-aaa5-4193eebcb79c
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=2, show_value=true)
end

# ╔═╡ 85db12e0-4146-49cd-b26a-736a1fc8467c
@bind c overlay_mask_bind(dilated_mask_3D)

# ╔═╡ f98b97fa-692d-477a-b6b4-7b72b262b684
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
        markersize=1.3,
        color=:red,
    )
    return fig
end

# ╔═╡ f6c5aab9-a8a2-485e-9d3c-cac9e47cee19
overlay_mask_plot(gaussian_poisson_calcification_3D_low_measure, dilated_mask_3D, c, "Overlayed Mask Low Energy")

# ╔═╡ 3776f9d2-41ea-4fa3-928c-12e0b85b3d90
overlay_mask_plot(gaussian_poisson_calcification_3D_high_measure, dilated_mask_3D, c, "Overlayed Mask High Energy")

# ╔═╡ Cell order:
# ╟─0a38e736-a860-49b6-9152-a51f33012117
# ╟─805556cc-0ed5-4b60-909b-933626c93e95
# ╠═f3363153-9287-482e-8049-bbc317aa8b51
# ╠═a77307a4-7d8e-4e40-8757-952a34d4d423
# ╟─d74d15b9-b330-4f31-af47-88a7f8d314b1
# ╠═080b868a-b8c2-445c-b53f-d587ac6e3288
# ╠═04b354df-36cd-422c-b5a1-09f19a0662ed
# ╠═dd93a5d6-fa5f-4a80-9b76-79abe755e833
# ╟─2f33b7b1-1a52-41e3-8521-f15f5594e971
# ╟─b8977e39-f5c4-4d3d-a6e2-104d5f11b549
# ╟─9cfd3199-e014-4c38-97d4-a2493e0ac437
# ╟─a528497a-0b2d-4561-836f-69a57d94e0db
# ╠═48ca6e64-010f-4460-8a2f-2eaf06cde1ec
# ╠═e27cf9d6-62b1-49d6-9bc1-07f7f6711af0
# ╟─8dbf0852-8139-45d1-acff-246202396715
# ╠═f5ca6a29-bbbb-4659-ab2e-ea93a5dfdafc
# ╠═239a4e84-48a8-48b2-a45b-32ed9d512ad5
# ╠═9cfbd6b5-3db1-489d-b57f-e7e7a48a0068
# ╟─25faf325-2ccb-47b2-80f8-a44882d53d8a
# ╠═dc58098b-8756-4119-9974-9335b5aeb5da
# ╠═cd417d4c-2c08-47d8-9925-f0eeaf063f63
# ╟─b6cdb1e7-be55-4898-91d2-624e6901965e
# ╠═a0cfcbd7-0f49-472d-a6e2-2de7b87bfb5e
# ╟─85db12e0-4146-49cd-b26a-736a1fc8467c
# ╟─f6c5aab9-a8a2-485e-9d3c-cac9e47cee19
# ╟─3776f9d2-41ea-4fa3-928c-12e0b85b3d90
# ╟─5a720fc1-d082-4f0b-b7ea-fbcc6a978345
# ╠═964b6d60-2e1b-468f-9d5c-7ce2c631e8c7
# ╠═68f05736-722c-49b8-9222-8c5548599cbf
# ╠═8b0dd679-6266-4ce7-960b-db56ddef1afa
# ╟─4384ae7f-ed4a-47ec-9979-1a7a25409fe7
# ╠═b4269f46-bdb4-4da6-acce-2f31f143d36a
# ╟─87c3f233-3607-4805-9394-16b6857a5887
# ╟─5d92ebf6-ef70-42d6-ab4b-3ce5d516f4fe
# ╟─296b5f5f-e8db-45a9-8010-3b98056d6781
# ╠═11e5c831-a033-4a7c-8019-8daebf5e92c3
# ╟─25c56212-e9e4-46fc-8cdd-f8d80fc7661d
# ╠═99e4b65f-c8bf-4a60-89d0-2f3d036a05dd
# ╠═744b47ef-2d5f-4f7b-989f-465d9ac7fc40
# ╠═0fb61ccb-43b7-4b1e-8b92-9692d0a42548
# ╟─eb39e77a-1256-4006-a625-b9dc4bd114eb
# ╠═6c22fc99-923a-486f-9d34-8c6d04f5e969
# ╟─6e652ad4-de3c-47cd-b8c6-9b0c2ebacee3
# ╠═1dd08bcb-e945-4e9c-9bbe-3fc37a843894
# ╠═3ab8f65a-ae68-49b0-b717-32de09d4b56f
# ╟─4c0bc41a-20bf-497d-8dbc-3b12389cb1aa
# ╟─321f8a2a-224e-48e2-9fb8-c8a8922c3e4a
# ╟─02721b3d-3a23-4633-9706-40c94856b91d
# ╟─c7e6d66d-18f9-4795-89ed-11c5c312d8f4
# ╟─c97c4a72-7857-42f7-abe4-8ad247ff538c
# ╟─b3d6ec94-4cc2-43d4-b64f-f4d6bf12d54d
# ╟─ee56071a-0d3d-48c3-9aff-73933c0c0146
# ╟─167385d2-8d62-4193-885e-4bbc3683461d
# ╟─6d6076eb-e1bc-44f3-a913-cda403af0d70
# ╟─4f031d7f-70c1-42d9-af64-8057178b7fca
# ╟─52acce68-1f6e-426b-ab84-c52e60b6870b
# ╟─66d311a5-0dd9-4d4e-99fd-377c3bfab82f
# ╟─cc9000bb-dc0c-4cdd-b15a-f2af91da5ea8
# ╟─f9c17967-6f97-44bb-9c8f-dbc35f09b0e2
# ╠═0a7c1547-beba-4562-9cd7-c2f882612182
# ╠═f394d548-3792-4999-bb54-5208e46e6a34
# ╠═b3a78f50-d709-40b8-a5f4-840534939091
# ╟─fa77e270-ab83-4e2f-a304-9f2ca7e44c0c
# ╠═a323a639-d40f-4f68-8e2b-822e6d70177c
# ╟─aaa53ab7-2aca-46c5-b15a-4506ec061ebd
# ╠═0a8f0a44-33fd-470f-b5fb-ca0138b21b5c
# ╠═d039d82e-3f26-4434-a9f8-cabad86b8771
# ╟─a81771c1-78ed-40b6-8517-ec71c5cb17a7
# ╟─6e33a03b-6484-49de-833c-7337b2cb1e2d
# ╟─939ef2cc-11b2-40ae-b489-59d15fe5a437
# ╟─de246ed9-912d-4c55-bd71-09972abb3d02
# ╟─653b3f92-e48e-4440-8913-fca7e90ca3fc
# ╟─2106147b-e6c6-48fd-8e65-ce1b329c06ae
# ╟─91b745c3-689f-4164-b3e5-4f49eecb93a4
# ╟─4a64c4c4-f713-4eb1-a18d-809e9816d594
# ╟─f026dcbf-0e74-49a0-94c0-32cfd77b869c
# ╟─bc224643-f914-42ee-96ea-ce5f40472350
# ╟─1add336d-dbcf-4466-a25a-043bde36a021
# ╟─4d8d6448-b60c-4072-959f-d32310c5a9dd
# ╟─a1b37552-936a-4d56-ba16-d280b29adaa8
# ╟─f95221cd-ae8d-4cf0-bfe9-cb7e005b242c
# ╟─40ef54ea-3afd-487b-8bea-806e33134b68
# ╟─147b5c20-b93d-45d8-aaa5-4193eebcb79c
# ╟─f98b97fa-692d-477a-b6b4-7b72b262b684

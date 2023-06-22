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
		Pkg.add("Statistics")
		Pkg.add("DICOM")
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
end

# ╔═╡ d6254d97-3b3b-4b2f-b7b3-8bc58e6f34c2
md"""
# Overview
This package contains four different state of the art coronary artery calcium (CAC) quantification methods.
1. Agatston scoring [1, 2]
2. Volume fraction calcium mass [3]
3. Integrated calcium mass [3, 4]
4. Spatially weighted calcium scoring [5, 6]

This notebook will briefly showcase how to use CalciumScoring.jl to quantify calcium using the gold standard Agatston scoring technique.

#### References
1. [Quantification of coronary artery calcium using ultrafast computed tomography](https://doi.org/10.1016/0735-1097(90)90282-t)
2. [Ultra-low-dose coronary artery calcium scoring using novel scoring thresholds for low tube voltage protocols—a pilot study ](https://doi.org/10.1093/ehjci/jey019)
3. [Coronary artery calcium mass measurement based on integrated intensity and volume fraction techniques](https://doi.org/10.1101/2023.01.12.23284482)
4. [Integrated intensity-based technique for coronary artery calcium mass measurement: A phantom study](https://doi.org/10.1002/mp.16326)
5. [An alternative method for quantifying coronary artery calcification: the multi-ethnic study of atherosclerosis (MESA)](https://doi.org/10.1186/1471-2342-12-14)
6. [Spatially Weighted Coronary Artery Calcium Score and Coronary Heart Disease Events in the Multi-Ethnic Study of Atherosclerosis](https://doi.org/10.1161/CIRCIMAGING.120.011981)

"""

# ╔═╡ f3a886e7-31e8-43dd-adc7-3fb50c7d7d23
md"""
## Import Packages
First, let's import the most up-to-date version of CalciumScoring.jl, which can be found on the main/master branch of the [GitHub repository](https://github.com/Dale-Black/CalciumScoring.jl). Be aware this can take a long time, especially if this is the first time being downloaded. Future work on this package will focus on improving this.
"""

# ╔═╡ df5551ce-01cc-49a2-b040-e1a2b6d794a4
TableOfContents()

# ╔═╡ fac54db6-e2a1-4956-8b7c-091938e74b8e
md"""
## Load Simulated Phantoms
Next, we need to load simulated images of cardiac phantoms.
"""

# ╔═╡ 29b04989-caa7-4e73-a9a1-cdad9e3509d8
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
end;

# ╔═╡ 0188bed1-93f4-4e95-b581-0ea64847faac
md"""
## Visualize Simulated Phantom
We loaded a simulated phantom containing various density calcifications. Let's visualize this image below along with a mask of the large, high-density, calcification.

This phantom is supposed to simulate a human thorax containing calcified coronary arteries. Since the phantom has identical geometry along each slice, the only difference along each slice will be from "quantum noise".
"""

# ╔═╡ 5c64de6c-2143-43a5-94cf-a7cdb6e4e210
@bind z1 PlutoUI.Slider(axes(dcm_array, 3); default = 2, show_value = true)

# ╔═╡ 0a3c60e5-b695-41e4-b43d-45b2b56c5b7e
let
	idxs = getindex.(findall(isone, mask_large_insert_high_density_3D[:, :, z1]), [1 2])
	idxs = [t[i] for t in idxs, i in 1:length(idxs[1])]
	
	f = Figure(resolution = (1200, 800))
	ax = Axis(
		f[1, 1],
		title = "DICOM Array"
	)
	heatmap!(transpose(dcm_array[:, :, z1]); colormap = :grays)
	scatter!(idxs[:, 2], idxs[:, 1]; markersize = 4, color = (:red, 0.2))
	hidedecorations!(ax)
	f
end

# ╔═╡ 6bb10787-21f4-48fc-bec2-e7c23ccf7a8b
md"""
# Agatston (Volume) Scoring
"""

# ╔═╡ 87286995-8337-410d-bc12-a2c1c638aff5
md"""
## Calculate Ground Truth Values
The goal with calcium scoring is to quantify the amount of calcium within a cardiac CT scan. Before we calculate the calcium via various calcium scoring techniques, lets calculate the known calcium mass. Since we know the density of the calcifications and the geometry of the inserts, we can calculate the ground truth calcium mass of each of the 9 calcium inserts. 

The calculation for volume is simple since the inserts are cylindrical
```math
\text{volume} = \pi \times r^2 \times \text{slice thickness} \times \text{number of slices}
\tag{1}
```

Mass is then simple to obtain from the volume and known density
```math
\text{mass} = \text{volume} \times \rho
\tag{2}
```

"""

# ╔═╡ f08c6b1e-0ae8-4b79-b282-d03597200a00
begin
	slice_thickness = dcms[1].meta[tag"Slice Thickness"] # mm
	pixel_spacing = dcms[1][tag"Pixel Spacing"] # mm
	voxel_spacing = [pixel_spacing..., slice_thickness]

	num_slices = length(axes(dcm_array, 3))
	ρ_high_density = 0.073 # mg/mm^3
	radius = (5/2) # mm
	calcium_volume_large = π * radius^2 * slice_thickness * num_slices
	calcium_mass_large = calcium_volume_large * ρ_high_density
end

# ╔═╡ 31ba16d9-e51f-48ef-8731-fce208de1c30
md"""
## Visualize Dilated Mask
Now let's dilate the mask a little bit, to account for partial volume effect (blurring) and visualize our ROI
"""

# ╔═╡ 83246419-b0e2-478c-bb0d-a7ab12c6ecd7
@bind z2 PlutoUI.Slider(axes(dcm_array, 3); show_value = true, default = 2)

# ╔═╡ 310ad8ab-838a-4080-ad09-076a56767174
md"""
## Results
"""

# ╔═╡ e48e6ae3-1fc4-4fa4-8d8c-52db5dd0c19d
md"""
# Agaston Mass Scoring
"""

# ╔═╡ 8b13e40c-b2fd-4f2d-856e-7605c5f5c0b7
md"""
We can also compute the calcium mass score, via the Agatston technique by computing the mass calibration factor. This factor translates an Agatston score to a mass score ([1](https://pubs.rsna.org/doi/epdf/10.1148/radiol.2432050808)).

We can assume the CT number of water is 0, and compute the correction factor by dividing the density of known calcium by the CT number of that same calcium. Let's load our simulated calibration phantom to acquire this.
"""

# ╔═╡ b6b0e32a-8823-413d-b3f9-16334ff715cd
md"""
## Load Calibration Rod
"""

# ╔═╡ eb718d32-2966-4869-8550-5b98ae8dc04a
begin
	dcm_path_calibration = joinpath(pwd(), "dcms/cal/200/30/80");
	dcms_calibration = dcmdir_parse(dcm_path_calibration)
	dcm_array_calibration = cat(
		[dcms_calibration[i][(0x7fe0, 0x0010)] for i in 1:length(dcms_calibration)]...;
		dims=3
	)
end;

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

# ╔═╡ c3837d06-7c00-49ca-8592-6b9e8e7ef528
center_insert1, center_insert2 = 187, 318

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
	hu_calcium, ρ_calcium = mean(calibration_rod), 0.200 # mg/mm^3
end

# ╔═╡ 2d3b75aa-da53-4490-a743-9a6ca405759f
mass_calibration = ρ_calcium / hu_calcium

# ╔═╡ 52c83cf7-8d82-4c84-ba57-3384ddec26f3
md"""
## Results
"""

# ╔═╡ cc1ad7c9-bfdf-43c3-bbf4-4e82112eb5af
md"""
# Next Steps

The core function, `score()`, implements all four types of CAC quantification methods. Agatston scoring is the current gold standard, but recent advances in the field of medical physics have introduced better approaches (in terms of accuracy and sensitivity) to calcium mass quantification. 

Check out the [Volume Fraction Calcium Mass](\volume-fraction.jl) tutorial to see how to implement one of these approaches.
"""

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
function dilate_recursively(mask, n)
    dilated_mask = copy(mask)
    for _ in 1:n
        dilated_mask = dilate(dilated_mask)
    end
    return dilated_mask
end

# ╔═╡ 51649204-97fe-460b-935e-06a29bb85d5b
begin
	dilated_mask_large_high_density = dilate_recursively(mask_large_insert_high_density, 2)
	dilated_mask_large_high_density_3D = cat(
		dilated_mask_large_high_density, 
		dilated_mask_large_high_density,
		dilated_mask_large_high_density;
		dims=3)
end;

# ╔═╡ 41340d94-d8e6-48cd-879f-707017969911
pure_calcification_large_high_density = create_mask(dcm_array, dilated_mask_large_high_density_3D);

# ╔═╡ 7e4c3269-2a2b-4b15-a373-4daf815a72ef
agatston_score, volume_score = score(pure_calcification_large_high_density, voxel_spacing, Agatston())

# ╔═╡ ade72d79-4082-4c4e-b14a-744651731d9d
_, _, mass_score = score(pure_calcification_large_high_density, voxel_spacing, mass_calibration, Agatston())

# ╔═╡ 320405f5-4feb-4c61-98f4-ff13643490ae
md"""
Let's compare that with the ground truth mass. We see that the ground truth mass = $(round.(calcium_mass_large, digits = 2)) is close to the calculated mass = $(round.(mass_score, digits = 2))
"""

# ╔═╡ 907fe9d2-a67e-4716-9c73-0ce429ae2b4f
let
	idxs = getindex.(findall(isone, dilated_mask_large_high_density_3D[:, :, z2]), [1 2])
	
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
	heatmap!(transpose(pure_calcification_large_high_density[:, :, z2]), colormap = :grays)
	f
end

# ╔═╡ Cell order:
# ╟─d6254d97-3b3b-4b2f-b7b3-8bc58e6f34c2
# ╟─f3a886e7-31e8-43dd-adc7-3fb50c7d7d23
# ╠═9419f973-985f-41af-b6a4-bbd6af9ea2a4
# ╠═df5551ce-01cc-49a2-b040-e1a2b6d794a4
# ╟─fac54db6-e2a1-4956-8b7c-091938e74b8e
# ╠═29b04989-caa7-4e73-a9a1-cdad9e3509d8
# ╟─0188bed1-93f4-4e95-b581-0ea64847faac
# ╟─5c64de6c-2143-43a5-94cf-a7cdb6e4e210
# ╟─0a3c60e5-b695-41e4-b43d-45b2b56c5b7e
# ╟─6bb10787-21f4-48fc-bec2-e7c23ccf7a8b
# ╟─87286995-8337-410d-bc12-a2c1c638aff5
# ╠═f08c6b1e-0ae8-4b79-b282-d03597200a00
# ╟─31ba16d9-e51f-48ef-8731-fce208de1c30
# ╠═51649204-97fe-460b-935e-06a29bb85d5b
# ╠═41340d94-d8e6-48cd-879f-707017969911
# ╟─83246419-b0e2-478c-bb0d-a7ab12c6ecd7
# ╟─907fe9d2-a67e-4716-9c73-0ce429ae2b4f
# ╟─310ad8ab-838a-4080-ad09-076a56767174
# ╠═7e4c3269-2a2b-4b15-a373-4daf815a72ef
# ╟─e48e6ae3-1fc4-4fa4-8d8c-52db5dd0c19d
# ╟─8b13e40c-b2fd-4f2d-856e-7605c5f5c0b7
# ╟─b6b0e32a-8823-413d-b3f9-16334ff715cd
# ╠═eb718d32-2966-4869-8550-5b98ae8dc04a
# ╟─f00fb7d4-5183-483a-90a2-73ff9be3a24b
# ╟─efd98ea7-b86a-41e1-9eca-6675dffc7edd
# ╠═c3837d06-7c00-49ca-8592-6b9e8e7ef528
# ╠═1870556b-2e19-4e64-9df1-338daae29805
# ╠═2d3b75aa-da53-4490-a743-9a6ca405759f
# ╟─52c83cf7-8d82-4c84-ba57-3384ddec26f3
# ╠═ade72d79-4082-4c4e-b14a-744651731d9d
# ╟─320405f5-4feb-4c61-98f4-ff13643490ae
# ╟─cc1ad7c9-bfdf-43c3-bbf4-4e82112eb5af
# ╟─982837e3-6d4c-4012-a014-886d38d19121
# ╟─107e9699-6baa-4978-a68f-9fc22d76d6a1
# ╟─ac9da8e2-aee9-4c79-8e7d-f35eb7e674f2
# ╟─4a14052c-5327-43a4-bf2e-ad739bf5213d
# ╟─140d7064-df9e-4dc8-886e-edb381792164
# ╟─78d99536-b253-491e-b08d-df7cb97a2063

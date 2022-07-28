### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 982c6570-6879-11eb-2af7-e1cd658b545f
begin
	using NeXLSpectrum
	using Gadfly
	md"""
# Fitting the Continuum
A workbook demonstrating how to fit a continuum model to a spectrum.
	"""
end

# ╔═╡ eff192d0-6879-11eb-2bd7-ed3157983e71
spec=loadspectrum(joinpath(@__DIR__,"..","test","ADM6005a spectra","ADM-6005a_1.msa"))

# ╔═╡ e00b5530-687a-11eb-0c4f-3145b609802e
begin
	set_default_plot_size(7inch,4inch)
	plot(spec, klms=[n"O",n"Al",n"Si",n"Ca",n"Ti",n"Zn",n"Ge"], xmax=12.0e3)
end

# ╔═╡ 21191060-6882-11eb-2a4e-d3fc8817c521
md"""
We need to know the spectrum composition to model the continuum.
"""

# ╔═╡ 774fb2b0-687b-11eb-0b06-9b3da85dc5c2
spec[:Composition]=parse(Material,"0.3398*O+0.0664*Al+0.0405*Si+0.0683*Ca+0.0713*Ti+0.1055*Zn+0.3037*Ge",name="ADM-6005a")

# ╔═╡ 0576efd0-6882-11eb-0252-c37d47ce32c4
md"""
Construct a detector that matches this spectrum with a resolution of 132.0 eV at 
Mn Kα and a low-level discriminator of 110 channels. The discriminator eliminates
the negative energy channels which contain the zero-strobe peak.
"""

# ╔═╡ 87d6a7b0-687b-11eb-0330-37d38343dd80
det=matching(spec, 132.0, 110)

# ╔═╡ 5cf17510-6881-11eb-3855-4fe7f874bc21
md"""
To fit the continuum, the model needs to know how the detector will respond to an 
X-ray. The detector response function (implemented as a matrix) is a function that convolves data collected on a perfect detector (perfect resolution and unity efficiency) to account for the resolution and the detector efficiency.
"""

# ╔═╡ 458d65a0-687c-11eb-18aa-29a72e160f55
resp = detectorresponse(det,SDDEfficiency(ModeledWindow(MoxtekAP33())))

# ╔═╡ 26bf4360-6885-11eb-316f-618f54f568b9
md"""
Fit a continuum model to the spectrum collected on `det` with response function `resp`.
"""

# ╔═╡ 3b2ffd10-6882-11eb-1445-19268242d6c2
brem = fittedcontinuum(spec,det,resp)

# ╔═╡ 14d20140-687d-11eb-21e4-8332aabd18b3
plot(spec,brem,klms=[n"O",n"Al",n"Si",n"Ca",n"Ti",n"Zn",n"Ge"],yscale=0.05, xmax=12.0e3)

# ╔═╡ 424745e0-687d-11eb-1668-cb12e7aa954e
plot(spec,brem,klms=[n"C",n"O",n"Al",n"Si",n"Ca",n"Ti",n"Zn",n"Ge"],yscale=0.05, xmax=2.0e3)

# ╔═╡ 16597580-6886-11eb-3675-b5dddf4ed005
md"""
The default Bremsstrahlung model is `Castellano2004a`.  You can specify alternative
models using the `brem` argument.  Available alternatives include `Kramers1923`, `Lifshin1974`, `Reed1975`, `Smith1975`, `Small1987`, `Trincavelli1997`, `Castellano2004a` and `Castellano2004b`.
"""

# ╔═╡ 774dab50-6885-11eb-36bd-df8a2a530275
brem2 = fittedcontinuum(spec,det,resp, brem=Kramers1923)

# ╔═╡ afced940-6885-11eb-01fc-114700a7753b
plot(spec,brem2,klms=[n"C",n"O",n"Al",n"Si",n"Ca",n"Ti",n"Zn",n"Ge"],yscale=0.05, xmax=12.0e3)

# ╔═╡ Cell order:
# ╟─982c6570-6879-11eb-2af7-e1cd658b545f
# ╠═eff192d0-6879-11eb-2bd7-ed3157983e71
# ╟─e00b5530-687a-11eb-0c4f-3145b609802e
# ╟─21191060-6882-11eb-2a4e-d3fc8817c521
# ╠═774fb2b0-687b-11eb-0b06-9b3da85dc5c2
# ╟─0576efd0-6882-11eb-0252-c37d47ce32c4
# ╠═87d6a7b0-687b-11eb-0330-37d38343dd80
# ╟─5cf17510-6881-11eb-3855-4fe7f874bc21
# ╠═458d65a0-687c-11eb-18aa-29a72e160f55
# ╟─26bf4360-6885-11eb-316f-618f54f568b9
# ╠═3b2ffd10-6882-11eb-1445-19268242d6c2
# ╠═14d20140-687d-11eb-21e4-8332aabd18b3
# ╠═424745e0-687d-11eb-1668-cb12e7aa954e
# ╟─16597580-6886-11eb-3675-b5dddf4ed005
# ╠═774dab50-6885-11eb-36bd-df8a2a530275
# ╠═afced940-6885-11eb-01fc-114700a7753b

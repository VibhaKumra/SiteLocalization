module SiteLocalization
using Images
using ImageFeatures
using ImageTransformations
using CoordinateTransformations
using StaticArrays
using Statistics
using Interpolations

export SNR, translateim, circle_add!, fit_circles, sigma_norm, rl_deconv, circle_inbounds, masked, stats, intp_polar

t = 2.262; # t-test p=0.05
function stats(distance, n_theta)
    n = length(distance)
    #avgdistance = sum(distance.*n_theta)./sum(n_theta)
    avgdistance = mean(distance)
    #stdev = std(distance, weights(n_theta), corrected=false)
    stdev = std(distance,corrected=false)
    t0 = avgdistance/(stdev./(n).^0.5)
    conf_int = t*stdev/n^0.5
    return avgdistance, stdev, n, t0, conf_int
end

"""
	circle_add!(img::AbstractArray{T, 2}, center, radius)

Adapted from ImageDraw (ellipse2d.jl). Draws a circle defined by `center` and `radius` to matrix `img`.
"""
function circle_add!(img::AbstractArray{T, 2}, center, radius) where T<:Integer
    ys = Int[]
    xs = Int[]
    cy, cx = center.I
    for i in cy : cy + radius
        for j in cx : cx + radius
            val = ((i - cy) / radius) ^ 2 + ((j - cx) / radius) ^ 2
            if val < 1
                push!(ys, i)
                push!(xs, j)
            end
        end
    end
    for (yi, xi) in zip(ys, xs)
        img[yi, xi] += 1
        if (yi != cy)
            img[2 * cy - yi, xi] += 1
        end
        if (xi != cx)
            img[yi, 2 * cx - xi] += 1
            if (yi != cy)
                img[2 * cy - yi, 2 * cx - xi] += 1
            end
        end
    end
    img
end

"""
	fit_circles(img, highp, lowp, blur,votingthres, mindistance, radius_range)

Fit circles using canny edge detection and circle hough transform.
"""
function fit_circles(img, highp, lowp, blur,votingthres, mindistance, radius_range)
    #display(typeof(img))
    #img = img .> 0.5
    img_edge = canny(img,  (Percentile(highp), Percentile(lowp)), blur)
    #println(extrema(img))
    #img_edge = img.> 0.01
    dx ,dy =imgradients(img, KernelFactors.ando5)
    img_phase = phase(dx,dy)
    centers, radii = hough_circle_gradient(img_edge, img_phase, radius_range, scale= 1, vote_threshold=votingthres,min_dist=mindistance)
    return centers, radii, img_edge
end

"""
	masked(img, center)

Find circle on original image, polar transformed coordinates
"""
function masked(img, center)
    transformation = PolarFromCartesian()
    mask = ((transformation(SVector(((i, j) .- center.I)...)), img[i,j]) for i=1:size(img, 1), j=1:size(img, 2))
    return mask
end


"""
	polar_transform(masked, radius)

Discretized polar transform.
"""
function polar_transform(masked, radius)
    circ = ceil(Int, 2*π*radius)
    polar_img = fill(-1., radius, circ)
    for (coord, val)=Iterators.filter(x -> x[1].r < radius, masked)
        i = 1+floor(Int, coord.r)
        j = 1+floor(Int, (π+coord.θ)*radius)
        polar_img[i, j] = val
    end
    return polar_img
end

"""
	duplicate_fill(polar_img)

Duplicate and fill in empty elements of polartransform
"""
function duplicate_fill(polar_img)
    ind = size(polar_img)
    polar_full = zeros(ind)
    circ = ind[2]
    for i = 1:ind[1]
        for j = 1:ind[2]
            nxt = 0
            val = -1.
            while val == -1.
                val = polar_img[i,1 + ((j+nxt+circ) % circ)]
                nxt = ~(nxt-signbit(nxt))
            end
            polar_full[i,j] = val
        end
    end
    return polar_full
end

"""
	rl_deconv(image::AbstractArray, psf::AbstractArray, iterations::Int)

Deconvolve `image` using the Richardson-Lucy deconvolution algorithm
"""
function rl_deconv(image::AbstractArray, psf::AbstractArray, iterations::Int)
    latent_est = Float64.(image)
    psf_hat = reflect(psf)
    interm = zeros(Float64, size(image))
    for i = 1:iterations
        latent_est.+= eps(Float64)
        # imfilter!(ArrayFireLibs(Algorithm.Mixed()), interm, latent_est, psf_hat)
        imfilter!(interm, latent_est, psf_hat)
        rel_blur = image ./ interm
        # imfilter!(ArrayFireLibs(Algorithm.Mixed()), interm, rel_blur, psf)
        imfilter!(interm, rel_blur, psf)
        latent_est .*= interm
    end
    return latent_est
end

"""
	circle_inbounds(center, radius, sz)

Test if the circle defined by `center` and `radius` is in-bounds for an array with size `sz`.
"""
function circle_inbounds(center, radius, sz)
    for i=1:length(sz)
        if (center[i]-radius < 1) || (center[i]+radius > sz[i])
            return false
        end
    end
    return true
end

function SNR(image)
    signal = maximum(image)
    noise = median(image)
    snr = 10*log10(signal/noise)
    return snr
end

abstranslation(x) = (-floor(Int, min(x, 0)),ceil(Int, max(x, 0)))
function translateim(image, d)
    pads = abstranslation.(d)
    warp(padarray(image, Pad(:replicate, (pads[1][1], pads[2][1]), (pads[1][2], pads[2][2]))), Translation(d))[indices_spatial(image)...]
end

function sigma_clipped_stats(img, sigma=3.0)
    nchanged = 1
    m = 0
    s = 0
    ori = sort(img[:])
    signal = ori
    bot, top = 1, length(ori)
    while nchanged != 0
        len = length(signal)
        mid = fld(len,2)
        m = if iseven(len)
            (signal[mid]+signal[mid+1])/2
        else
            signal[mid+1]
        end

        s = std(signal)
        nbot, ntop = searchsortedfirst(ori, m-sigma*s), searchsortedlast(ori, m+sigma*s)
        nchanged = abs(bot-nbot)+abs(top-ntop)
        bot, top = nbot, ntop
        sig = @view ori[bot:top]
        signal = sig
    end
    return m, s
end

function sigma_norm(img)
    img = Float64.(img)
    m, s = sigma_clipped_stats(img)
    img .-= m + 3.0*s
    img = img ./ maximum(img)
    clamp01.(img)
end

function read_avg_from_txt(txtfile)
    avg_text = map(row -> split(row, ',')[2], readlines(txtfile)[1:end-2])
    return parse.(Float64, map(text -> split(text, ':')[2], avg_text))
end

function datastats(data, filtfun)
    println(length(data))
    data = filter(filtfun, data)
    err = std(data)/sqrt(length(data))
    println(mean(data))
    println(err)
    return(data, err)
end

datastats(data) = datastats(data, x->!isnan(x))


"""
    intp_polar(img, radius, angles, center; intp = BSpline(Constant()))

Perform polar transform with interpolation of `img` centered in `center`.
"""
function intp_polar(img, radius, angles, center; intp = BSpline(Constant()))
    interpolator = interpolate(img, intp)
    [interpolator(center[1]+r*sin(theta),center[2]+r*cos(theta))::Float64 for r in 0.0:1:radius-1.0, theta in angles]
end
    

end

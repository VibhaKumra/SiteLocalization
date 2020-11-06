# -*- coding: utf-8 -*-
using Images
using ImageFeatures
using CoordinateTransformations
using StaticArrays

#FUNCTIONS
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

# Adopted from ImageDraw (ellipse2d.jl)
# Add circle to matrix
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

function highpass(img, f)
    fftimg = fftshift(fft(ifftshift(gray.(img)), 1:2))
    sx, sy = size(fftimg)
    for x = 1:sx, y=1:sy
        if (x - sx/2.0)^2 + (y - sx/2.0)^2 < f^2
            fftimg[x,y] = 0
        end
    end
    return abs.(fftshift(ifft(ifftshift(fftimg))))
end


#Function for fitting circles using canny edge detector
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

#Function for finidng circle on original image, polar transformed coordinates
function masked(img, center)
    transformation = PolarFromCartesian()
    mask = ((transformation(SVector(((i, j) .- center.I)...)), img[i,j]) for i=1:size(img, 1), j=1:size(img, 2))
    return mask
end

# Discretized polartransform
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

#Function for duplicating and filling in empty elements of polartransform
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

# Function for Richardson Lucy Deconvolution
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

# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---


include("functionsSiteLoc.jl")

using Images
using ImageDraw
using Plots
using StatsPlots
using Distributions
using ImageFiltering
using StatsBase
using DSP
using Glob
using ProgressMeter
using RegisterQD
using Photometry
using ImageTransformations
using Distances
using Interpolations
using Primes
using ImageIO


image = reshape(1:100, (10,10))
abstranslation(x) = (-floor(Int, min(x, 0)),ceil(Int, max(x, 0)))
function translateim(image, d)
    pads = abstranslation.(d)
    warp(padarray(image, Pad(:replicate, (pads[1][1], pads[2][1]), (pads[1][2], pads[2][2]))), Translation(d))[indices_spatial(image)...]
end
translateim(image, (-1.6, -2.8))

padarray(image, Pad(:replicate, (1,2), (3,4)))

# +
function sigma_clipped_stats(img, sigma=3.0)
    nchanged = 1
    m = 0
    s = 0
    ori = sort(img[:])
    signal = oriantb_shift,
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




gr() #set plotting backend
Localisationdata= "testruns"


open("RESULTS_$Localisationdata.txt", "w") do f

    #BioformatsLoader.init(memory= 2048)
    path = "artifacts"

    #Loading Images
    image_dirs = [""] # add image directory
    get_imgs = [fname for d in image_dirs for fname in glob("*.nd2", d)] #read files
    
    # Set Variables
    pixel_size = 30.8/1.5# in nanometers
    shiftwindow = 3 #for sliding average of peak alignment, 0 if only aligning max of peaks
    maxpeakswindow = 0 # sliding average for sum of several circles
    margin = 15 #how much to add on radius of fitted circle
    rmin = 24 #radius min in pixels
    rmax= 30 #radius max
    radius_range = rmin:rmax
    padded_size(radius) = 2*radius+2
    
    #For canny edge detection and circle hough transform
    highp = 99  #in percentile
    lowp = 95  #in percentile
    blur = 3.5 #blur for segmentation
    votingthres = 24 #in pixels, sensitivity of circle detection
    mindistance = 50 #in pixels, the minimum distance between two circle centers
    
    # For deconvolution
    PSF = 4 # measured or known sigma of psf
    psf = Kernel.gaussian([PSF,PSF],[21,21])
    psf = psf./maximum(psf)
    iterations = 0 #no of iterations of deconvolution algorithm
    
    #ALIGNMENT
    al_shift = (0.4, -0.15) # adjust for measured chromatic aberration
    
    #Setting various variables
    max_radius = (1 + maximum(radius_range) + margin)*2 #set max size for plots
    println("max_radius = $max_radius")
    allmemb = zeros(padded_size(max_radius))
    allantb = zeros(padded_size(max_radius))

    distance = Array{Array{Float64,1},1}()
    n_theta = Array{Array{Int,1},1}()
    PerBac_distance = Array{Array{Float64,1},1}()
    PerBac_n_theta = Array{Array{Int,1},1}()

    # Find max peak for alignment using a sliding average
    window_mean(window,v) = [mean(v[i-window:i+window]) for i=1+window:length(v)-window]
    movmax(window,v) = window + argmax(window_mean(window, v))
    intp = BSpline(Cubic(Free(OnGrid())))
    intp = BSpline(Constant())
    
    # Function for polar transform 
    function intp_polar(img, radius, angles, center)
        interpolator = interpolate(img, intp)
        [interpolator(center[1]+r*sin(theta),center[2]+r*cos(theta))::Float64 for r in 0.0:1:radius-1.0, theta in angles]
    end
    
    #Following code run for all imgs in loaded folder
    for	file_index = 4 #1:Int(length(get_imgs))
        fname = (get_imgs[file_index])
        if !isfile(fname)
            continue
        end
        
        data_fname = "$path/artifacts/$fname.txt"
        channels, _ = load(query(fname))
        firstSNRantb = SNR(channels[1][2,:,:])
        firstSNRmemb = SNR(channels[1][1,:,:])
        if firstSNRantb <= 0.03 || firstSNRmemb <= 0.03
           continue
        end
        println(firstSNRantb)
        println(firstSNRmemb)
        push!(distance, Array{Float64,1}())
        push!(n_theta, Array{Int,1}())
        Timeframes = 0
        @showprogress 1 fname for t = 1
            antibodies = translateim(channels[t][2,:,:], al_shift) |> sigma_norm
            membrane = channels[t][1,:,:] |> sigma_norm
            
            currentSNRantb = SNR(antibodies)
            currentSNRmemb = SNR(membrane)

            if currentSNRantb <= 0.3*firstSNRantb || currentSNRmemb <= 0.3*firstSNRmemb # SNR threshold for time series
                break
            end
            
            Timeframes = t
            #Deconvolve images
            membrane = rl_deconv(membrane, psf, iterations)
            antibodies = rl_deconv(antibodies, psf, iterations)

            #Normalizing images
            my_norm(image) = (e -> @. (image-e[1]+eps())*(1/(e[2]+eps()-e[1])))(extrema(image))
            antibodies = my_norm(antibodies)
            membrane = my_norm(membrane)
            bkgantb= median(antibodies)
            stdbkgantb= std(antibodies)

            #implementing fit_circles on membrane channel
            centers, radii, membraneEdge = fit_circles(membrane, highp, lowp, blur, votingthres, mindistance, radius_range)
            # check if circle is in bound
            circles = Iterators.filter(x -> circle_inbounds(x[1], x[2]+margin, size(membrane)), zip(centers, radii))
            centers2, radii2, antibodyEdge = fit_circles(antibodies,highp, lowp, blur, votingthres, mindistance, radius_range)

            #drawing circles on membrane Image
            membraneEdge = Gray.(membraneEdge)
            for (center, radius) in circles
                draw!(membraneEdge, CirclePointRadius(center, radius))
            end

            #Save Images
            save("$path/membrane_edge.png", membraneEdge)
            save("$path/membrane_norm.png", membrane)
    

            # circle fill: on zero matrix += ones for each circle, to localize overlap of circles
            circle_img = zeros(Int,size(membrane))
            for (center, radius) in circles
                try
                    circle_add!(circle_img, center, radius+margin)
                finally
                end
            end

            # Find all elements in circle_img > 1
            elements_overlap = findall(x -> x > 1, circle_img)
            list_overlap = CartesianIndices(circle_img)
            list_overlap = list_overlap[elements_overlap]

            # Find elements in circle
            in_circle(c, r, x, y) = ((x-c[1])^2 + (y-c[2])^2) <= r^2
            totalmemb = zeros(padded_size(max_radius))
            totalantb = zeros(padded_size(max_radius))
            sumdistance = zeros(1)
            n = zeros(1)
            
            #Extract mask in both images
            for (center, radius) in circles
                radius += margin

                memb_masked = masked(membrane, center)
                antb_masked = masked(antibodies, center)

                transformation = PolarFromCartesian()
                circ = ceil(Int, 2*π*radius)
                overlap = Iterators.filter(t -> in_circle(center, radius, t[1], t[2]), list_overlap)
                polar_overlap = (transformation(SVector((Tuple(t - center)))) for t=overlap)
                discrete_theta = Set(1+floor(Int, (π+c.θ)*radius) for c=polar_overlap)

                crop(center) = center-radius:center+radius
                memb_polar_full = clamp01.(intp_polar(membrane[crop(center[1]), crop(center[2])], radius, 0:1/radius:2*pi, (radius+1, radius+1)))
                antb_polar_full = clamp01.(intp_polar(antibodies[crop(center[1]), crop(center[2])], radius, 0:1/radius:2*pi, (radius+1, radius+1)))
                save("$path/memb_polar_full.png", memb_polar_full)
                save("$path/antb_polar_full.png", antb_polar_full)
                replace!(memb_polar_full, NaN =>0)
                replace!(antb_polar_full, NaN =>0)


                radius = size(memb_polar_full, 1)
                n_radius = zeros(padded_size(radius))
                r = 1:radius
                
                maxPeaks = mapslices(x -> movmax(shiftwindow, x) , memb_polar_full, dims=(1,))[:]
                designmethod = Butterworth(5)
                ff = digitalfilter(Lowpass(0.1),designmethod)
           
                maxPeaks = clamp.(round.(Int,filtfilt(ff, maxPeaks)), 1, radius)

                memb_shift = zeros(padded_size(radius), circ)
                antb_shift = zeros(padded_size(radius),circ)

                for θ = 1:circ
                    if !in(θ, discrete_theta)
                        
                        shift_ind = @. (1:radius) + radius + 1 - maxPeaks[θ]
                        memb_shift[shift_ind,θ] = memb_polar_full[:,θ]
                        antb_shift[shift_ind,θ] = antb_polar_full[:,θ]
                        n_radius[shift_ind] .+= 1
                    end
                end
                replace!(memb_shift, NaN =>0)
                replace!(antb_shift, NaN =>0)
                save("$path/memb_shift.png", memb_shift)
                save("$path/antb_shift.png", antb_shift)
                
                antb_shift = Iterators.filter(x -> mean(x) > 0, eachcol(antb_shift))
                try
                    antb_shift = sum(antb_shift)
                catch ArgumentError
                    continue
                end
                memb_shift = mapslices(sum, memb_shift, dims=(2,))

                range_total = @. (-radius:radius+1) + 1 + max_radius
                totalantb[range_total] += antb_shift./n_radius
                totalmemb[range_total] += memb_shift./n_radius

                rad = collect(-radius:radius)

                #Sliding average
                maxPeaksmemb = mapslices(x -> movmax(maxpeakswindow, x), memb_shift, dims=(1,))
                

                #Sliding average
              
                #findlocalmaxima
                maxPeaksantb = map(x -> movmax(maxpeakswindow, x), [antb_shift[memb_peak:end,i] for (i, memb_peak) in enumerate(maxPeaksmemb)])
                function findrightmaxima(x)
                    maximas = [c.I[1] for c in findlocalmaxima(x)]
                    m = mean(x)*antb_shift
                    first(c for c in reverse(maximas) if x[c] >= m)
                end
                maxPeaksantb = maxPeaksmemb .+ maxPeaksantb
                memb_shift, antb_shift = my_norm(memb_shift), my_norm(antb_shift)
                p = plot([ (1:size(memb_shift,1))*pixel_size],
                    [memb_shift], xlims=[(1+0.5*radius)*pixel_size,(1+1.5*radius)*pixel_size])
                display(p)
                savefig(p, "$path/Intensityplotex_WGAsumaligned.svg")
              
                append!(distance[end], maxPeaksantb - maxPeaksmemb)
                append!(n_theta[end], circ - length(discrete_theta))

            end #(center, radius) in circle
            
            allantb .+= totalantb
            allmemb .+= totalmemb
            replace!(allmemb, NaN =>0)
            
            replace!(allantb, NaN =>0)
        end #t
        avgdistance, stdev, n, t0, conf_int = stats(distance[end], n_theta[end])
        if n < 2 && false
            break
        end
        display("avg:$(avgdistance.*pixel_size) nm,confidence interval:$(conf_int.*pixel_size) nm, std:$(stdev.*pixel_size) nm, n:$n, t0:$(t0), dpts = $avgdistance")
        write(f, "Stats:, avg:$(avgdistance*pixel_size), std:$(stdev*pixel_size), n:$(n), t0:$(t0), conf_int:$(conf_int*pixel_size)\n")
    end #file_index

    totalmean = mean(mean.(distance))
    totalstdev = std(mean.(distance))/sqrt(Int(length(get_imgs)))
    println(totalmean*pixel_size)
    println(totalstdev*pixel_size)
    #Save  txt file of parameters, stats and results
    write(f, "totalmean: $(totalmean*pixel_size), totalstderror: $(totalstdev*pixel_size)\n")
    write(f, "Parameters:,  pixel_size:$pixel_size, shiftwindow:$shiftwindow,  maxpeakswindow:$maxpeakswindow, margin:$margin, radius_range:$radius_range, highp:$highp, lowp:$lowp, blur:$blur, votingthres:$votingthres, mindistance:$mindistance\n")
    write(f, "$distance")
    
end


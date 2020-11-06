import Pkg
Pkg.activate(".")
using SiteLocalization
using Images
using ImageDraw
using Plots
using StatsPlots
using Distributions
using ImageFiltering
using CoordinateTransformations
using StatsBase
using DSP
using Glob
using ProgressMeter
using RegisterQD
using Distances
using FileIO
import GR
GR.inline("svg")
gr() #set plotting backend


path = "artifacts"

# Loading Images
if length(ARGS) > 0
    get_imgs = ARGS
else
    get_imgs = [""] # add images
end

# Set Variables
pixel_size = 30.8/1.5# in nanometers
shiftwindow = 7 #for sliding average of peak alignment, 0 if only aligning max of peaks
maxpeakswindow = 0 # sliding average for sum of several circles
margin = 20 #how much to add on radius of fitted circle
rmin = 24 #radius min in pixels
rmax= 30 #radius max
radius_range = rmin:rmax
padded_size(radius) = 2*radius+2

# For canny edge detection and circle hough transform
highp = 98.8  #in percentile
lowp = 99.9  #in percentile
blur = 3.5 #blur for segmentation
votingthres = 23 #in pixels, sensitivity of circle detection
mindistance = 50 #in pixels, the minimum distance between two circle centers

# For deconvolution
psf_sigma = 4 # measured or known sigma of psf
psf = Kernel.gaussian([psf_sigma,psf_sigma],[21,21])
psf = psf./maximum(psf)
iterations = 0 #no of iterations of deconvolution algorithm

# ALIGNMENT
al_shift = (0.4, -0.15) # adjust for measured chromatic aberration

# Setting various variables
max_radius = (1 + maximum(radius_range) + margin)*2 #set max size for plots

localisation_data = "testruns"
open("RESULTS_$localisation_data.txt", "w") do f

    allmemb = zeros(padded_size(max_radius))
    allantb = zeros(padded_size(max_radius))

    distance = Array{Array{Float64,1},1}()
    n_theta = Array{Array{Int,1},1}()
    PerBac_distance = Array{Array{Float64,1},1}()
    PerBac_n_theta = Array{Array{Int,1},1}()

    # Find max peak for alignment using a sliding average
    window_mean(window,v) = [mean(v[i-window:i+window]) for i=1+window:length(v)-window]
    movmax(window,v) = window + argmax(window_mean(window, v))

    #Following code run for all imgs in loaded folder
    for	file_index = 1:Int(length(get_imgs))
        fname = get_imgs[file_index]
        
        data_fname = "$path/artifacts/$fname.txt"
        channels = [load(file) for file in glob(fname)]
        if length(channels) == 0
            println("Found no matching files for: $fname! Skipping...")
        end
        if length(channels) == 1 && ndims(channels[1]) != 3
            channels = channels[1]
        end
        firstSNRantb = SNR(channels[1][2,:,:])
        firstSNRmemb = SNR(channels[1][1,:,:])
        if firstSNRantb <= 0.03 || firstSNRmemb <= 0.03
           continue
        end
        push!(distance, Array{Float64,1}())
        push!(n_theta, Array{Int,1}())
        @showprogress 1 fname for t = 1:length(channels)
            antibodies = translateim(channels[t][2,:,:], al_shift) |> sigma_norm
            membrane = channels[t][1,:,:] |> sigma_norm
            
            currentSNRantb = SNR(antibodies)
            currentSNRmemb = SNR(membrane)

            if currentSNRantb <= 0.3*firstSNRantb || currentSNRmemb <= 0.3*firstSNRmemb # SNR threshold for time series
                break
            end
            
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
        println("avg:$(avgdistance.*pixel_size) nm,confidence interval:$(conf_int.*pixel_size) nm, std:$(stdev.*pixel_size) nm, n:$n, t0:$(t0), dpts = $avgdistance")
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

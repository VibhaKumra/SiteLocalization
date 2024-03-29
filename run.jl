import Pkg
Pkg.activate(".")
using Dates
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
using StaticArrays
using DelimitedFiles
using JSON3
using ArgParse
import Profile
import GR
GR.inline("svg")
gr() #set plotting backend


path = "artifacts"

function get_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--pixelsize"
            help = "image pixel size in nanometers"
            arg_type = Float64
            default = 30.8/1.5
        "--shiftwindow"
            help = "for sliding average of peak alignment, 0 if only aligning max of peaks"
            arg_type = Int
            default = 7
        "--maxpeakswindow"
            help = "sliding average for sum of several circles"
            arg_type = Int
            default = 0
        "--iterations"
            help = "RL iterations"
            arg_type = Int
            default = 0
        "--psf_sigma"
            help = "Measured sigma of psf, for Richardsson Lucy deconvolution"
            arg_type = Int
            default = 4
        "--highp"
            help = "For canny edge detection and circle hough transform, in percentile"
            arg_type = Float64
            default = 99.8
        "--lowp"
            help = "For canny edge detection and circle hough transform, in percentile"
            arg_type = Float64
            default = 95.0
        "--blur"
            help = "blur for segmentation"
            arg_type = Float64
            default = 3.5
        "--mindistance"
            help = "in pixels, the minimum distance between two circle centers"
            arg_type = Int
            default = 50
        "--margin"
            help = "fitted circle radius margin"
            arg_type = Int
            default = 20
        "--votingthreshold"
            help = "sensitivity of circle detection in pixels"
            arg_type = Int
            default = 20
        "--rmin"
            help = "fitted circle minimum radius"
            arg_type = Int
            default = 36
        "--rmax"
            help = "fitted circle maximum radius"
            arg_type = Int
            default = 44
        "--pxmult"
            help = "pixel interpolation"
            arg_type = Float64
            default = 2.0
        "--al-shifty"
            help = "adjust for measured chromatic aberration"
            arg_type = Float64
            default = -0.63
        "--al-shiftx"
            help = "adjust for measured chromatic aberration"
            arg_type = Float64
            default = -1.35
        "--first-snr-limit"
            help = "SNR limit for first frame for analyzing a time sequence"
            arg_type = Float64
            default = 0.0
        "--snr-ratio-limit"
            help = "SNR limit in terms of ratio to SNR of first frame at which to stop analyzing a time sequence"
            arg_type = Float64
            default = 0.3
        "--limit-outside"
            help = "limit search for peaks to outside of reference channel"
            action = :store_true
        "--maxdist"
            help = "maximum distance from reference channel within which to search for peaks"
            arg_type = Int
            default = 0
        "--debug-plots"
            help = "enable saving of plots and images for debugging purposes"
            action = :store_true
        "--enable-sigma-norm"
            help = "enable saving of plots and images for debugging purposes"
            action = :store_true
        "--name"
        "images"
            help = "images to analyze"
            required = true
            nargs = '+'
            arg_type = String
    end
    parse_args(s)
end

padded_size(radius) = 2*radius+2
function main(f, csvf, get_imgs; debug_plots, pxmult, limit_outside, pixel_size,
        shiftwindow, maxpeakswindow, margin, radius_range, highp, lowp, path,
        enable_sigma_norm, blur, votingthres, mindistance, psf, iterations,
        al_shift, max_radius, maxdist, snr_ratio_limit, first_SNR_limit, rmin)
    
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
            println("Found no matching files for: $(fname)! Skipping...")
        end
        if length(channels) == 1 && ndims(channels[1]) != 3
            channels = channels[1]
        end
        channels = [channels[1][:,:,2i-1:2i] for i in 1:fld(size(channels[1],3), 2)]
        firstSNRantb = SNR(channels[1][:,:,2])
        firstSNRmemb = SNR(channels[1][:,:,1])
        
        if firstSNRantb <= first_SNR_limit || firstSNRmemb <= first_SNR_limit
           continue
        end
        push!(distance, Array{Float64,1}())
        push!(n_theta, Array{Int,1}())
        @showprogress 1 fname for t = 1:length(channels)
            antibodies, membrane = let frame = channels[t]
                ab = translateim(frame[:,:,2], al_shift)
                memb = frame[:,:,1]
                if enable_sigma_norm
                    sigma_norm(ab), sigma_norm(memb)
                else
                    ab, memb
                end
            end
            
                
            currentSNRantb = SNR(antibodies)
            currentSNRmemb = SNR(membrane)
            
            if currentSNRantb < snr_ratio_limit*firstSNRantb || currentSNRmemb < snr_ratio_limit*firstSNRmemb
                break
            end
            
            if debug_plots && t == 1
                save("$path/antb_norm/$(basename(fname)).png", antibodies)
            end
                
            #Deconvolve images
            membrane = rl_deconv(membrane, psf, iterations)
            membrane ./= maximum(membrane)
            antibodies = rl_deconv(antibodies, psf, iterations)
            antibodies ./= maximum(antibodies)
                
            if debug_plots && t == 1
                save("$path/membrane/$(basename(fname)).png", membrane)
            end

            #Normalizing images
            my_norm(image) = (e -> @. (image-e[1]+eps())*(1/(e[2]+eps()-e[1])))(extrema(image))
            antibodies = my_norm(antibodies)
            membrane = my_norm(membrane)
            bkgantb= median(antibodies)
            stdbkgantb= std(antibodies)
                
            
            padw = margin
            antibodies = padarray(antibodies, Pad(:replicate, (padw, padw), (padw, padw))).parent
            membrane = padarray(membrane, Pad(:replicate, (padw, padw), (padw, padw))).parent

            #implementing fit_circles on membrane channel
            centers, radii, membraneEdge = fit_circles(membrane, highp, lowp, blur, votingthres, mindistance, radius_range)
                
            # check if circle is in bound
            circles = Iterators.filter(x -> circle_inbounds(x[1], x[2]+margin, size(membrane)), zip(centers, radii))

            #drawing circles on membrane Image
            membraneEdge = Gray.(membraneEdge)
            for (center, radius) in circles
                draw!(membraneEdge, CirclePointRadius(center, radius))
            end

            #Save Images
            if debug_plots && t == 1
                save("$path/membrane_edge/$(basename(fname)).png", membraneEdge)
                save("$path/membrane_norm/$(basename(fname)).png", membrane)
                save("$path/antb_norm_deconv/$(basename(fname)).png", antibodies)
            end
    

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
            for (center, circle_radius) in circles
                radius = circle_radius + margin

                memb_masked = masked(membrane, center)
                antb_masked = masked(antibodies, center)

                transformation = PolarFromCartesian()
                circ = ceil(Int, 2*π*radius)
                overlap = Iterators.filter(t -> in_circle(center, radius, t[1], t[2]), list_overlap)
                polar_overlap = (transformation(SVector((Tuple(t - center)))) for t=overlap)
                discrete_theta = Set(1+floor(Int, (π+c.θ)*radius) for c=polar_overlap)

                crop(center) = center-radius:center+radius
                function perform_interpolation(channel)
                    cropped = channel[crop(center[1]), crop(center[2])]
                    clamp01nan.(intp_polar(cropped, radius, 0:1/radius:2*pi, (radius+1, radius+1); pxmult=pxmult))
                end
                
                memb_polar_full = perform_interpolation(membrane)
                antb_polar_full = perform_interpolation(antibodies)
                    
                if debug_plots && t == 1
                    save("$path/memb_polar_full/$(basename(fname)).png", memb_polar_full)
                    save("$path/antb_polar_full/$(basename(fname)).png", antb_polar_full)
                end


                radius = size(memb_polar_full, 1)
                n_radius = zeros(padded_size(radius))
                r = 1:radius
                
                maxPeaks = mapslices(x -> movmax(shiftwindow, x) , memb_polar_full, dims=(1,))[:]
                #[floor(Int,(rmin)*pxmult):end,:]
                #designmethod = Butterworth(5)
                #ff = digitalfilter(Lowpass(0.1),designmethod)
           
                #maxPeaks = clamp.(round.(Int,filtfilt(ff, maxPeaks)), 1, radius)

                memb_shift = zeros(padded_size(radius), circ)
                antb_shift = zeros(padded_size(radius), circ)

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
                if debug_plots && t == 1
                    save("$path/memb_shift/$(basename(fname)).png", memb_shift)
                    save("$path/antb_shift/$(basename(fname)).png", antb_shift)
                end
                
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
                    
                maxPeaksantb = if limit_outside
                    maxPeaksantb = map(x -> movmax(maxpeakswindow, x) - 1, [antb_shift[memb_peak-offset:end,i] for (i, memb_peak) in enumerate(maxPeaksmemb)])
                    maxPeaksmemb .+ maxPeaksantb
                elseif maxdist > 0
                    range(memb_peak) = max(1,floor(Int, memb_peak-maxdist*pxmult)):min(size(antb_shift,1), ceil(Int, memb_peak+maxdist*pxmult))
                    map(enumerate(maxPeaksmemb)) do (i, memb_peak)
                        x = antb_shift[range(memb_peak),i]
                        movmax(maxpeakswindow, x) - 1 + minimum(range(memb_peak))
                    end
                else
                    #Sliding average
                    mapslices(x -> movmax(maxpeakswindow, x), antb_shift, dims=(1,))
                end
                
                

                function findrightmaxima(x)
                    maximas = [c.I[1] for c in findlocalmaxima(x)]
                    m = mean(x)*antb_shift
                    first(c for c in reverse(maximas) if x[c] >= m)
                end
                    
                memb_shift, antb_shift = my_norm(memb_shift), my_norm(antb_shift)

                if debug_plots && t == 1
                    p = plot(
                            [ (1:size(memb_shift,1))*pixel_size, (1:size(memb_shift,1))*pixel_size ],
                            [memb_shift, antb_shift],
                            #xlims=[(1+0.5*radius)*pixel_size,(1+1.5*radius)*pixel_size],
                            xlims=[(0)*pixel_size,1700],
                            labels=["Membrane" "Antibody"],
                            legend=false,
                            lw=3,
                            tickfontsize =12
                    )

                    #vline!([maxPeaksmemb, maxPeaksantb]*pixel_size, labels=["Membrane" "Antibody"])
                    #annotate!([((maxPeaksantb+maxPeaksmemb)/2*pixel_size, 0.2, "$(maxPeaksantb-maxPeaksmemb) px")])

                    #w_img = plot(p, plot(Gray.(imadjustintensity(membrane))),plot(Gray.(imadjustintensity(antibodies))), size=(1500,500), layout=(1,3))
                    mkpath("$path/Intensityplotex_WGAsumaligned/")
                    
                    savefig(p, "$path/Intensityplotex_WGAsumaligned/$(basename(fname)).svg")
                end
              
                append!(distance[end], (maxPeaksantb - maxPeaksmemb))
                append!(n_theta[end], circ - length(discrete_theta))
                let crops = (crop(center[1]), crop(center[2]))
                    writedlm(f, [fname t distance[end][end] n_theta[end][end] "$center" circle_radius minimum(crops[1]) maximum(crops[1]) minimum(crops[2]) maximum(crops[2])])
                    flush(f)
                end

            end #(center, radius) in circle
            
            allantb .+= totalantb
            allmemb .+= totalmemb
            replace!(allmemb, NaN =>0.0)
            replace!(allantb, NaN =>0.0)
        end #t
        avgdistance, stdev, n, t0, conf_int = stats(distance[end], n_theta[end])
        if n < 2 && false
            break
        end
        println("avg:$(avgdistance.*pixel_size) nm,confidence interval:$(conf_int.*pixel_size) nm, std:$(stdev.*pixel_size) nm, n:$n, t0:$(t0), dpts = $avgdistance")
        write(f, "Stats:, avg:$(avgdistance*pixel_size), std:$(stdev*pixel_size), n:$(n), t0:$(t0), conf_int:$(conf_int*pixel_size)\n")
        writedlm(csvf, [fname avgdistance stdev avgdistance*pixel_size stdev*pixel_size n t0 conf_int*pixel_size])
        flush(csvf)
    end #file_index
    mean_values = filter(x-> !isnan(x) && x > -100/pixel_size && x< 150/pixel_size, mean.(distance))
    display(mean_values)
    totalmean = mean(mean_values)
    totalstdev = std(mean_values)
    println(totalmean*pixel_size)
    println(totalstdev*pixel_size)
    #Save  txt file of parameters, stats and results
    write(f, "totalmean: $(totalmean*pixel_size), totalstderror: $(totalstdev*pixel_size)\n")
    write(f, "Parameters:,  pixel_size:$pixel_size, shiftwindow:$shiftwindow,  maxpeakswindow:$maxpeakswindow, margin:$margin, radius_range:$radius_range, highp:$highp, lowp:$lowp, blur:$blur, votingthres:$votingthres, mindistance:$mindistance\n")
    write(f, "$distance")
end

function run_with_args(args; name="")
    get_imgs = args["images"]
    offset=24;
    # Set Variables
    pixel_size = args["pixelsize"]/args["pxmult"] # in nanometers
    shiftwindow = args["shiftwindow"] #for sliding average of peak alignment, 0 if only aligning max of peaks
    maxpeakswindow = args["maxpeakswindow"] # sliding average for sum of several circles
    margin = args["margin"] #how much to add on radius of fitted circle
    rmin = args["rmin"] #radius min in pixels
    rmax = args["rmax"] #radius max
    enable_sigma_norm = args["enable-sigma-norm"] #radius max
    radius_range = rmin:rmax

    # For canny edge detection and circle hough transform
    highp = args["highp"]  #in percentile
    lowp = args["lowp"]  #in percentile
    blur = args["blur"] #blur for segmentation
    votingthres = args["votingthreshold"] #in pixels, sensitivity of circle detection
    mindistance = args["mindistance"] #in pixels, the minimum distance between two circle centers

    # For deconvolution
    psf_sigma = args["psf_sigma"] # measured or known sigma of psf
    psf = Kernel.gaussian([psf_sigma,psf_sigma],[21,21])
    psf = psf./maximum(psf)
    iterations = args["iterations"] #no of iterations of deconvolution algorithm

    # ALIGNMENT
    #al_shift = (0.4, -0.15) # adjust for measured chromatic aberration
   
    #al_shift = (-0.63, -1.35) # adjust for measured chromatic aberration
    al_shift = args["al-shifty"], args["al-shiftx"]

    # Setting various variables
    max_radius = ceil(Int, (1 + maximum(radius_range) + margin)*args["pxmult"]) #set max size for plots

    debug_plots = args["debug-plots"]
    pxmult = args["pxmult"]
    limit_outside = args["limit-outside"]

    timestamp = "$(Dates.format(now(), dateformat"YY-mm-dd_HHMMSS"))"

    out_dir = "results"
    mkpath(out_dir)
    artpath = mkpath(joinpath(path, name))

    open(joinpath(out_dir, "$(timestamp)_$(name)_params.json"), "w") do io
        JSON3.write(io, args)
    end

    localisation_data = "testruns"
    open(joinpath(out_dir, "$(timestamp)_$(name)_results.tsv"), "w") do csvf
        writedlm(csvf, ["fname" "avgdistance" "stdev" "avgdistance_micron" "stdev_micron" "n" "t0" "conf_int_micron"])
        open("RESULTS_$localisation_data.txt", "w") do f
            main(f, csvf, get_imgs; debug_plots, pxmult, limit_outside, pixel_size, maxdist=args["maxdist"],
                shiftwindow, maxpeakswindow, margin, radius_range, highp, lowp, enable_sigma_norm,
                blur, votingthres, mindistance, psf, iterations, al_shift, max_radius,
                snr_ratio_limit=args["snr-ratio-limit"], first_SNR_limit=args["first-snr-limit"], rmin, path = artpath
               )
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()

    run_with_args(args)
end

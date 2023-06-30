using ChemStats, DataFrames, CSV, Plots
#PART 1 - After SAFD, CompCreate and Feature alignment

#Convert from string of vectors into vectors
function getVec(matStr) 

    if matStr[1] .== '[' && !contains(matStr,"\"")
        if contains(matStr, ", ")
            str = split(matStr[2:end-1],", ")
        else
            str = split(matStr[2:end-1]," ")
        end
    elseif contains(matStr,"\"")
        matStr = matStr[3:end-2]
        tempstr= split(matStr,"\"")
        str = []
        for i = 1:length(tempstr)
            if tempstr[i] != ", " 
                tempstr1 = split(tempstr[i][2:end-1],", ")
                for k = 1:length(tempstr1)
                    if tempstr1[k] != ""
                        push!(str,tempstr1[k])
                    else 
                        tempstr1[k] = "0"
                        push!(str,tempstr1[k])
                    end
                end
            end
        end

    elseif matStr[1] .== 'A'
        if contains(matStr, ", ")
            str = split(matStr[5:end-1],", ")
        else
            str = split(matStr[5:end-1]," ")
        end
    elseif matStr[1] .== 'F'
        if matStr .== "Float64[]"
            return []
        else
            str = split(matStr[9:end-1],", ")
        end
    else
        println("New vector start")
    end
    if (length(str) .== 1 && cmp(str[1],"") .== 0) 
        return []
    else
        str = parse.(Float64, str)
        return str
    end
end
#Extract fragments from feature alignment reports
function no_frags(data1)
    fragments= DataFrame()
    for i = 1:size(data1,1)

        str = data1[i]
        if str!="0"
            frag= getVec(str)
            numb = length(frag)
            append!(fragments,[DataFrame(Nr = i) DataFrame(Fragments = [frag]) DataFrame(N_frags = numb)])
        elseif str == "0"
            frag = Vector{Float64}([0])
            numb = Vector{Float64}([0])
            append!(fragments,[DataFrame(Nr = i) DataFrame(Fragments= [frag]) DataFrame(N_frags = numb)])
        end
    end
    return fragments
end
#Extract intensity from feature alignment reports
function frag_inten(data2)
    inten= DataFrame()
    for i = 1:size(data2,1)

        str = data2[i]
        if str!="0"
            
            int= getVec(str)
            append!(inten,[DataFrame(Nr = i) DataFrame(Intensity = [int])])
        elseif str == "0"
        int = Vector{Float64}([0])
        append!(inten,[DataFrame(Nr = i) DataFrame(Intensity = [int])])
        end
    end
    return inten
end
#Create one dataframe of each run
function change(rep)
    mz= rep.MS2Comp
    int = rep.MS2CompInt
    feats = zeros(1:length(rep.MeasMass))
    for r= 1:size(rep.MeasMass,1)
        values = rep.MeasMass[r]
        rounds= round(values, digits =3)
        feats[r] = rounds
    end

        df = innerjoin(no_frags(mz), frag_inten(int), on = :Nr)
        df[!,:Feature] = feats
    return df
end

pathin_mz = "C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\clean_mz" #files with m/z values
filename_mz = readdir(pathin_mz, join=true)

pathin_int = "C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\clean_int" # files with fragment intensity
filename_int = readdir(pathin_int, join = true)


rep02 = CSV.read("C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\clean_nonrep\\170426_02_Edith_120417_CCF_16_Cent_report_comp.CSV", DataFrame)
rep03= CSV.read("C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\clean_nonrep\\170426_03_Edith_120417_CCF_17_Cent_report_comp.CSV", DataFrame)
rep14= CSV.read("C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\clean_nonrep\\170425_14_Edith_120417_CCF_14_Cent_report_comp.CSV", DataFrame)


function changedf(num)
    df_mz= DataFrame()
    df_int = DataFrame()
    df_mz= CSV.read(filename_mz[num],DataFrame)
    df_int = CSV.read(filename_int[num],DataFrame)
    r_mz1= no_frags(df_mz[:,10])
    r_int1= frag_inten(df_int[:,10])
    r_mz2 = no_frags(df_mz[:,11])
    r_int2 = frag_inten(df_int[:,11])
    df = innerjoin(r_mz1, r_int1, r_mz2, r_int2, on = :Nr, makeunique=true)
    df.Features = df_mz.AveMass
    return df
end

rep1 = changedf(1)
rep10 = changedf(2)
rep11=  changedf(3)
rep12=  changedf(4)
rep13= changedf(5)
rep4 =  changedf(6)
rep5 =  changedf(7)
rep6 = changedf(8)
rep7= changedf(9)
rep8=  changedf(10)
rep9= changedf(11)

function sum_int(x)
    means_int = zeros(size(x,1))
    for j = 1:size(x,1)
        means_ints=mean(x[j,1:2])
        means_int[j] = means_ints
    end
    return means_int
end

function nonunique1(v)
    sort_vec = sort(v)
    uni = unique(@view sort_vec[[diff(sort_vec).==0; false]])
    return uni
end # gives repeating numbers in a vector

function data4(rep)
    combinespect = DataFrame()
   
    for k=1:size(rep,1)
        ms2_1= [rep.Fragments[k] rep.Intensity[k]]
        ms2_2= [rep.Fragments_1[k] rep.Intensity_1[k]]
        index_1= findall(in(ms2_2[:,1]), ms2_1[:,1])
        index_2= findall(in(ms2_1[:,1]), ms2_2[:,1])
        feat = round(rep.Features[k], digits=3)
    
        if length(index_1) > length(index_2)
            vect = rep.Fragments[k] 
            frag_rep = nonunique1(vect)# gives repeating fragment
            index_rep = findall(vect->vect==frag_rep[1], vect) # get index of the repeated fragments
            int_remove = minimum(ms2_1[1:end.∈ [index_rep],2]) # get the min intensity of the repeated fragment
            intens = ms2_1[:,2]
            index_rep_in = findall(intens->intens==int_remove,intens) # get index of the min intensity

            ms2_1 = ms2_1[Not([index_rep_in[1]]), :] # remove the repeated fragment of min intensity
            index_1= findall(in(ms2_2[:,1]), ms2_1[:,1]) # change index_1
        elseif length(index_2) > length(index_1)
            vec = rep.Fragments_1[k] 
            frag_rep = nonunique1(vec)
            index_rep = findall(vec->vec==frag_rep[1], vec)
            int_rep =maximum(ms2_2[1:end.∈ [index_rep],2])
            int_remove = ms2_2[1:end.∈ [index_rep],2]
            int_remove = minimum(ms2_2[1:end.∈ [index_rep],2])
            intens = ms2_2[:,2]
            index_rep_in = findall(intens->intens==int_remove,intens)

            ms2_2 = ms2_2[Not([index_rep_in[1]]), :]
            index_2= findall(in(ms2_1[:,1]), ms2_2[:,1])
        else 
            ms2_1= ms2_1
            ms2_2= ms2_2
        end

        if (index_1 == Int64[] || index_2 == Int64[]) || (ms2_1 == ms2_2 && ms2_1[1,1] == 0)
            mz_val= cat(ms2_1[:,1], ms2_2[:,1], dims=1) 
            mz_int = cat(ms2_1[:,2], ms2_2[:,2], dims=1) 
            mz_val = mz_val[mz_val .!= 0]
            mz_int= mz_int[mz_int .!= 0]
            n_frag = length(mz_val)
            append!(combinespect,[DataFrame(FragmentsMZ = [mz_val]) DataFrame(FragmentsInt = [mz_int]) DataFrame(N_frags = [n_frag]) DataFrame(Feature= [feat])])
            
        elseif rep.Fragments[k] == rep.Fragments_1[k]
            same_int = [ms2_1[1:end .∈ [index_1],2] ms2_2[1:end .∈ [index_2],2]]# gives same mz values diff intensity
            mz_val= ms2_1[1:end .∈ [index_1],1]# gives same mz values  
            mz_int= sum_int(same_int)
            n_frag = length(mz_val)
            append!(combinespect,[DataFrame(FragmentsMZ = [mz_val]) DataFrame(FragmentsInt = [mz_int]) DataFrame(N_frags = [n_frag]) DataFrame(Feature= [feat])])
           
        else
            diff_mz = cat(ms2_1[1:end.∉ [index_1],1], ms2_2[1:end.∉ [index_2],1], dims=1)
            diff_int =cat(ms2_1[1:end.∉ [index_1],2], ms2_2[1:end.∉ [index_2],2], dims=1)
            same_int = [ms2_1[1:end .∈ [index_1],2] ms2_2[1:end .∈ [index_2],2]]# gives same mz values diff intensity
            same_mz= ms2_1[1:end .∈ [index_1],1]# gives same mz values  
            means_int = sum_int(same_int)
            same_mz_int = [same_mz means_int]
            mz_val= cat(same_mz, diff_mz, dims =1)
            mz_int= cat(means_int, diff_int, dims =1)
            n_frag = length(mz_val)
            append!(combinespect,[DataFrame(FragmentsMZ = [mz_val]) DataFrame(FragmentsInt = [mz_int]) DataFrame(N_frags = [n_frag]) DataFrame(Feature= [feat])])
            
        end
    end
    return combinespect
end

spec1 = data4(rep1)
spec4 = data4(rep4)
spec5 = data4(rep5)
spec6 = data4(rep6)
spec7= data4(rep7)
spec8 = data4(rep8)
spec9 = data4(rep9)
spec10 = data4(rep10)
spec11 = data4(rep11)
spec12 = data4(rep12)
spec13 = data4(rep13)
spec2 = change(rep02)
spec3 = change(rep03)
spec14 = change(rep14)


#Part 2 - Finding Common Features between Runs

#alignment between runs - removing empty cells and non common features
pathin_feat = "C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\spect align comp"
filename_feat= readdir(pathin_feat, join = true)

    data =CSV.read(filename_feat[1],DataFrame)
    columns = occursin.("Cent_report",names(data))
    alignedfeaut  = data[:,columns]
    indeces_not0 = []
    indeces_0 = []
    for k = 1:size(data,1)
        if alignedfeaut[k,1] .!= "0" && alignedfeaut[k,2] .!= "0" && alignedfeaut[k,3].!= "0" && alignedfeaut[k,4] .!= "0" 
            push!(indeces_not0,k)
        else
            push!(indeces_0,k)
        end
    end
    newdata = data[indeces_not0,:] 
    CSV.write("C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\spect align comp",newdata)

#recombining replicates - 
    function data1(rep)
    combinespect = DataFrame()
   
    for k=1:size(rep,1)
        ms2_1= rep.Fragments[k]
        ms2_2= rep.Fragments_1[k]
        index_1= findall(in(ms2_2), ms2_1)
        index_2= findall(in(ms2_1), ms2_2)
        feat = round(rep.Features[k], digits=3)

        if (index_1 == Int64[] || index_2 == Int64[]) || (ms2_1 == ms2_2 && ms2_1[1,1] == 0)
            mz_val= cat(ms2_1, ms2_2, dims=1) 
            mz_val = mz_val[mz_val .!= 0]
            n_frag = length(mz_val)
            append!(combinespect,[DataFrame(Feature= [feat]) DataFrame(FragmentsMZ = [mz_val])  DataFrame(N_frags = [n_frag])])
            
        elseif rep.Fragments[k] == rep.Fragments_1[k]
            same_int = [ms2_1[1:end .∈ [index_1],2] ms2_2[1:end .∈ [index_2],2]]# gives same mz values diff intensity
            mz_val= ms2_1[1:end .∈ [index_1],1]# gives same mz values  
            append!(combinespect,[DataFrame(Feature= [feat]) DataFrame(FragmentsMZ = [mz_val]) DataFrame(N_frags = [n_frag])])
           
        else
            diff_mz = cat(ms2_1[1:end.∉ [index_1],1], ms2_2[1:end.∉ [index_2],1], dims=1)
            same_mz= ms2_1[1:end .∈ [index_1],1]# gives same mz values  
            mz_val= cat(same_mz, diff_mz, dims =1)
            n_frag = length(mz_val)
            append!(combinespect,[DataFrame(Feature= [feat]) DataFrame(FragmentsMZ = [mz_val]) DataFrame(N_frags = [n_frag])])
            
        end
    end
    return combinespect
end

#For runs that don't have replicates
function different(data1,feats)
    fragments= DataFrame()
    for i = 1:size(data1,1)

        str = data1[i]
        if str!="0"
            frag= getVec(str)
            numb = length(frag)
            featur = feats[i]
            append!(fragments,[DataFrame(Feature= featur) DataFrame(Fragments = [frag]) DataFrame(N_frags = numb)])
        elseif str == "0"
            frag = Vector{Float64}([0])
            numb = Vector{Float64}([0])
            featur = feats[i]
            append!(fragments,[DataFrame(Feature= featur)  DataFrame(Fragments= [frag]) DataFrame(N_frags = numb)])
        end
    end
    return fragments

end

function replicates(reps,rep1,rep2)
    r_mz1= no_frags(rep1)
    r_mz2 = no_frags(rep2)
    df = innerjoin(r_mz1, r_mz2, on = :Nr, makeunique=true)
    df.Features = reps.AveMass
return df
end

#Getting Common Features and Fragments between Runs 
pathin_align = "C:\\Users\\Lisa\\OneDrive - UvA\\Minor data\\Serum data\\Clean spectra\\spect align comp\\clean_aligned"
filename_align= readdir(pathin_align, join = true)

#Comparing different runs according to accumulation time or swath windows

#acc 20, 40, 60 with 60windows
    win60 =CSV.read(filename_align[5],DataFrame)

      
    rep008 = replicates(win60, win60[:,11],win60[:,12])
    spec8 = data1(rep008)
    rename!(spec8,:FragmentsMZ => :Spec8,  :N_frags =>:N_frags8)

    spec2 =different(win60[:,11], win60.AveMass)
    spec2 =rename!(spec2,:Fragments => :Spec2,  :N_frags =>:N_frags2)
    spec3 = different(win60[:,12], win60.AveMass)
    spec3 =rename!(spec3,:Fragments => :Spec3, :N_frags =>:N_frags3)

    acums= [20 40 60]
    frags_win60= [sum(spec3.N_frags3)  sum(spec8.N_frags8) sum(spec2.N_frags2)]
   
    show(names(win80))
#acc 20, 40, 60 with 80windows
    win80 =CSV.read(filename_align[6],DataFrame)
    rep009 = replicates(win80, win80[:,11],win80[:,14])
    spec9= data1(rep009)
    rename!(spec9,:FragmentsMZ => :Spec9,  :N_frags =>:N_frags9)
    rep001= replicates(win80, win80[:,10],win80[:,13])
    spec1= data1(rep001)
    rename!(spec1,:FragmentsMZ => :Spec1,  :N_frags =>:N_frags1)
    spec14 =different(win80[:,12], win80.AveMass)
    spec14 =rename!(spec14,:Fragments => :Spec14,  :N_frags =>:N_frags14)

    acums= [20 40 60]
    frags_win80= [sum(spec14.N_frags14)  sum(spec1.N_frags1) sum(spec9.N_frags9)]
    
 #acc 20, 40, 60 with 100windows
    show(names(win100))
    win100 =CSV.read(filename_align[4],DataFrame)
    rep006 = replicates(win100, win100[:,10],win100[:,13])
    spec6= data1(rep006)
    rename!(spec6,:FragmentsMZ => :Spec6,  :N_frags =>:N_frags6)
    rep011= replicates(win100, win100[:,12],win100[:,15])
    spec11= data1(rep011)
    rename!(spec11,:FragmentsMZ => :Spec11,  :N_frags =>:N_frags11)
    rep007 = replicates(win100, win100[:,11],win100[:,14])
    spec7= data1(rep007)
    rename!(spec7,:FragmentsMZ => :Spec7,  :N_frags =>:N_frags7)
   
    acums= [20 40 60]
    frags_win100= [sum(spec6.N_frags6)  sum(spec11.N_frags11) sum(spec7.N_frags7)]
    
 #win 60, 80, 100 with 20ms
    acc20 =CSV.read(filename_align[1],DataFrame)
    show(names(acc20))
    rep006 = replicates(acc20, acc20[:,10],acc20[:,13])
    spec6= data1(rep006)
    rename!(spec6,:FragmentsMZ => :Spec6,  :N_frags =>:N_frags6)
    spec14 =different(acc20[:,11], acc20.AveMass)
    spec14 =rename!(spec14,:Fragments => :Spec14,  :N_frags =>:N_frags14)
    spec3 =different(acc20[:,12], acc20.AveMass)
    rename!(spec3,:Fragments => :Spec3,  :N_frags =>:N_frags3)

    wins= [60 80 100]
    frags_acc20= [sum(spec3.N_frags3)  sum(spec14.N_frags14) sum(spec6.N_frags6)]

    #win 60, 80, 100 with 40ms
    acc40 =CSV.read(filename_align[2],DataFrame)
    show(names(acc40))
    rep08 = replicates(acc40, acc40[:,11],acc40[:,14])
    spec8= data1(rep08)
    rename!(spec8,:FragmentsMZ => :Spec8,  :N_frags =>:N_frags8)
    
    rep01 = replicates(acc40, acc40[:,10],acc40[:,13])
    spec1= data1(rep01)
    rename!(spec1,:FragmentsMZ => :Spec1,  :N_frags =>:N_frags1)

    rep11 = replicates(acc40, acc40[:,12],acc40[:,15])
    spec11= data1(rep11)
    rename!(spec11,:FragmentsMZ => :Spec11,  :N_frags =>:N_frags11)

    wins= [60 80 100]
    frags_acc40= [sum(spec8.N_frags8)  sum(spec1.N_frags1) sum(spec11.N_frags11)]

    #win 60, 80, 100 with 60ms
    acc60 =CSV.read(filename_align[3],DataFrame)
    show(names(acc60))
    rep009 = replicates(acc60, acc60[:,11],acc60[:,14])
    spec9= data1(rep009)
    rename!(spec9,:FragmentsMZ => :Spec9,  :N_frags =>:N_frags9)
    rep007= replicates(acc60, acc60[:,10],acc60[:,13])
    spec7= data1(rep007)
    rename!(spec7,:FragmentsMZ => :Spec7,  :N_frags =>:N_frags7)
    spec2 =different(acc60[:,12], acc60.AveMass)
    rename!(spec2,:Fragments => :Spec2,  :N_frags =>:N_frags2)

    wins= [60 80 100]
    frags_acc60= [sum(spec2.N_frags2) sum(spec9.N_frags9)  sum(spec7.N_frags7) ]

    frags_acc60
#Bar plots
p1 = bar([(20,frags_win60[1]), (40,frags_win60[2]),(60,frags_win60[3])],xticks = 20:20:60, fill=["pink", "lightblue", "lightgreen"], legend = false, ylabel = "Total Number of Fragments", title="60")
p2 = bar([(20,frags_win80[1]), (40,frags_win80[2]),(60,frags_win80[3])], xticks = 20:20:60, fill=["pink", "lightblue", "lightgreen"], legend = false, title="80", ms=2, ma=0.5, xlabel = "Acumulation Time (ms)")
p3 = bar([(20,frags_win100[1]), (40,frags_win100[2]),(60,frags_win100[3])], xticks = 20:20:60,fill=["pink", "lightblue", "lightgreen"],legend = false, title="100")
display("image/png", plot(p1,p2,p3,layout = (1,3)))
xticks!(0:20:60)

pa = bar([(60,frags_acc20[1]), (80,frags_acc20[2]),(100,frags_acc20[3])],xticks = 60:20:100, fill=["pink", "lightblue", "lightgreen"], legend = false, ylabel = "Total Number of Fragments", title="20ms")
pb= bar([(60,frags_acc40[1]), (80,frags_acc40[2]),(100,frags_acc40[3])], xticks = 60:20:100, fill=["pink", "lightblue", "lightgreen"], legend = false, title="40ms", ms=2, ma=0.5, xlabel = "Number of SWATH Windows")
pc = bar([(60,frags_acc60[1]), (80,frags_acc60[2]),(100,frags_acc60[3])], xticks = 60:20:100,fill=["pink", "lightblue", "lightgreen"],legend = false, title="60ms")
display("image/png", plot(pa,pb,pc,layout = (1,3)))


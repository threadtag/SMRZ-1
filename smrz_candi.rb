# Author: Hengyi Jiang <hengyi.jiang@gmail.com>
# purpose: to pick the smrz ribozyme from the candidate list 
# based on smrz_candi_pick.rb
# usage: ruby smrz_candi.rb [options] input_file
# 2020-04-12
# last update 2021-04-26

require 'strscan'
require 'optparse'
require 'fileutils'
require_relative("fasta.rb")
require_relative('seq_utils.rb')
require_relative("rna_fold.rb")
require_relative('header_parser.rb')

options={}
warning_on=false
temperature=50
opt_parser = OptionParser.new do |opt|
    opt.banner = "Usage: ruby smrz_candi.rb [OPTIONS] fasta_file"

    opt.on("-l","--maxinner l",Numeric,"maximum inner loop length") do |l|
        options[:l] = l.to_i
    end

    opt.on("-m","--minimumI m",Numeric,"minimum stemI length") do |m|
        options[:m] = m.to_i
    end

    opt.on("-n","--minimumII n",Numeric,"minimum stemII length") do |n|
        options[:n] = n.to_i
    end

    opt.on("-s","--strict","strict mode") do 
        options[:s]=true
    end

    opt.on("-d","--drawing","default with drawing") do
        options[:d]=true
    end

    opt.on("-t","--temperature t","folding temperature") do |t|
        temperature=t.to_i
    end

    opt.on("-w","--warnings","with warnings") do |w|
        warning_on=true
    end

    opt.on("-c","--countonly","count") do |c|
        options[:c]=true
    end
    opt.on("-r","--repress","repress positive output") do 
        options[:r]=true
    end
    opt.on("-i","--info","report hit information") do 
        options[:i]=true
    end

end.parse!



unless options.has_key?(:m)
   options[:m]=4
end

unless options.has_key?(:n)
    options[:n]=2
end


def warning(warning_on, msg)
    if warning_on
        puts msg
    end
end



if (ARGV.size <1 )
    puts "usage:ruby smrz_candi.rb [OPTIONS] input_file"
    puts "input_file: the fasta file name of the testing candidate"
    puts "Options:"
    puts "-l: inner loop maximum length"
    puts "-m: minimum stemI length"
    puts "-n: minimum stemII length,default 2"
    puts "-d: output drawing"
    puts "-s: strict mode"
    puts "-w: show rejection warning for debugging"
    puts "-c: count only for preliminary candidate with paired ends, no further test was done"
    puts "-t: folding temperature"
    puts "-i: report hit type, for debugging"
    puts "-r: repress positive hit"
	abort
end

file = ARGV.shift
if options.has_key?(:l)
    pattern =/T[AGT]CTACG\w(\w{1,#{options[:l]})\w[ACG]CGGGG/
            # i  +1 23456 7                            54321 -
else
    pattern = /T[AGT]CTACG\w(\w{1,})\w[ACG]CGGGG/
end
#          $1    $2        $3      $4      $5   $6         $7  
           
if !FileTest.exist?(file)
    puts  "file #{file} not found"
    abort
end

starting_seqs=Fasta_collection.new(file)
starting_seqs.load!

filtered_seqs=Fasta_collection.new(nil)
filtered_mids=Fasta_collection.new(nil)
matching_seqs=Fasta_collection.new(nil)


stem_loop_seq_count=0
filtered_seqs_count=0

# stage 1 collect candidates with preliminary sequence matching
starting_seqs.each do |fasta|
    fasta.sanitize
    if m=fasta.seq.match(/#{pattern}/)
        idx = fasta.seq.index(m[0])
        idx_r =idx+m[0].size-1
    else
        warning(warning_on,"#{fasta.header}:rejected from primary pattern match")
        next
    end
    if idx<1
        warning(warning_on,"#{fasta.header}:rejected in idx aquiring ")
        next
    end
    seq_head = fasta.seq[0..(idx-1)]
    seq_tail = fasta.seq[(idx_r+1)..-1]

   # first test the stemI check
    stemI_len=rev_comp_length(seq_head,seq_tail,true)+1
    if stemI_len < options[:m]
        warning(warning_on,"#{fasta.header}:stemI len ")
        next
    end
    if options.has_key?(:c) and options[:c]==true
        stem_loop_seq_count +=1        
    end

    core_idx=idx
    core_len= m[0].size
    special_A_idx = core_idx +7            
    special_C_idx = core_idx+core_len-7
 

    special_A = fasta.seq[special_A_idx]
    special_C = fasta.seq[special_C_idx]
    mid_loop=m[1]


    # second test special pairs 
    unless ["AC","AG","AA","AT","GT","CC","CT","CG","TG"].include?(special_A+special_C)
        warning(warning_on,"#{fasta.header}:rejected for special pair")
        next
    end
    if options.has_key?(:c) and options[:c]==true
        filtered_seqs_count +=1        
    else
        filtered_seqs.add!(fasta)
        filtered_mids.add!(Fasta_seq.new(fasta.header,mid_loop))
    end
end
starting_seqs=nil # free memory
if options.has_key?(:c) and options[:c]==true
    puts "all possible step loop candidates: #{stem_loop_seq_count}, after special pair filtering: #{filtered_seqs_count}"
    abort
end


#stage 2: fold all the prilimary candidates, global fold and local fold
folded_rnas=Folded_collection.new
folded_rnas.fold_fasta_collection!(filtered_seqs,temperature)

folded_mids=Folded_collection.new
folded_mids.fold_fasta_collection!(filtered_mids,temperature)

# puts folded_rnas.size
# puts folded_mids.size

#starg 3 : pick candidate by taking reference to the folded structures
filtered_seqs.each do |fasta|
    m=fasta.seq.match(/#{pattern}/)
    idx = fasta.seq.index(m[0])
    idx_r =idx+m[0].size-1    
    
    drawing = ""
    # pattern =/T[AGT]CTACG\w(\w{1,#{options[:l]})\w[ACG]CGGGG/
    #         # i  +1 23456 7                            54321 -
    core_idx=idx
    core_len= m[0].size
    core_last_idx=core_idx+core_len-1
    tripleI_A_idx = core_idx+4
    tripleI_C_idx = core_idx+5
    tripleI_G_idx = core_idx+6
    special_A_idx = core_idx +7
            
    special_C_idx = core_idx+core_len-7
    tripleII_A_idx  = core_idx+core_len-6
    tripleII_C_idx  = core_idx+core_len-5
    tripleII_G_idx  = core_idx+core_len-4

    stemI_left_1=core_idx-1
    stemI_left_2=core_idx-2

    stemI_right_1=core_last_idx+1
    stemI_right_2=core_last_idx+2

    a7_idx=core_idx+1
    a24_idx=tripleII_A_idx

    special_A = fasta.seq[special_A_idx]
    special_C = fasta.seq[special_C_idx]
    a7=fasta.seq[a7_idx]
    a24=fasta.seq[tripleII_A_idx]
    seq_head = fasta.seq[0..(idx-1)]
    seq_tail = fasta.seq[(idx_r+1)..-1]

    mid_loop=m[1]
    stemI_len=rev_comp_length(seq_head,seq_tail)+1
   

    stemII_folding_ok = false
    fold_ok = false
    output_tag =""

    if special_A=="A" and special_C == "C"
        special_AC_perfect=true
    else
        special_AC_perfect =false
    end

    if a7=="A"
        a7_perfect = true
    else
        a7_perfect = false
    end

    if a24=="A"
        a24_perfect =true
    else
        a24_perfect=false
    end


    # first test globe folding
    rna_full = folded_rnas.get(fasta.header)
    if rna_full ==nil
        # p fasta
        next
    end

    unless ["C","G"].include?(fasta.seq[stemI_left_1])
        # this paring should be strong
        warning(warning_on,"#{fasta.header}:stemI C:G pairing")
        next
    end

    # exception list
    enabled_mid_list=["ACATGGACGA"]

    #rejecting some special cases
    if mid_loop.size >5 and mid_loop.size < 10
        # puts "#{rna_full.unpaired_between?(special_A_idx+2,special_C_idx-2)},#{fasta.seq[special_A_idx..(special_A_idx+1)]}"
        if a7=="T" and fasta.seq[special_A_idx..(special_A_idx+1)]=="CG" and fasta.seq[(special_C_idx-1)..special_C_idx]=="CG"
            #case INVZ-13
            warning(warning_on,"#{fasta.header}:rejected for alternative folding enforced by other CG pairs")
            next
        end
    end


    if a24=="G" and special_C=="C" and rna_full.paired?(special_A_idx+1,special_C_idx-1)
        #case INV-Z9, Z10 
        strech_left=rna_full.get_stretch(special_A_idx+1)        
       
        if strech_left !=nil and (strech_left[0]-strech_left[1]).abs >=4
           #pass
        else
            warning(warning_on,"rejection:#{fasta.header}, possible alternative folding")
            next
        end
    end
 

    # special case testing
    check_comp = max_self_complementary(mid_loop)
    rna_mid = folded_mids.get(fasta.header)

    if mid_loop.size > 5
       
        # test if it is cannonical case 
        if rna_full.paired?(tripleI_C_idx,tripleII_G_idx) and 
            rna_full.paired?(tripleI_G_idx,tripleII_C_idx)  and 
            rna_full.paired?(special_A_idx+1,special_C_idx-1) and
            rna_full.unpaired?([tripleI_A_idx,tripleII_A_idx,special_A_idx,special_C_idx]) 
            #connonical , Tr9 case
                     
            unless ["CG","GC"].include?(special_A+a24)
                # special_A and A24 pairing up block structure               
                fold_ok =true
                output_tag="#1"                
            end            
            
        elsif a24=="G" and rna_full.unpaired?([tripleI_A_idx,special_A_idx]) and rna_full.paired?(special_A_idx+1,special_C_idx-1)
            unless rna_full.unpaired?([tripleI_C_idx,tripleI_G_idx])
                if rna_full.paired?(tripleI_C_idx,tripleII_A_idx) and rna_full.paired?(tripleI_G_idx,special_C_idx)
                    fold_ok =true
                    output_tag="#2"
                end
            end
                        
        
        elsif ["T","G"].include?(a7) and rna_full.unpaired?([tripleI_A_idx,tripleII_A_idx,special_A_idx]) and  rna_full.paired?(special_A_idx+1,special_C_idx-1)
            if rna_full.paired?(tripleI_G_idx,special_C_idx) and rna_full.unpaired?([tripleII_C_idx,tripleII_C_idx])
                fold_ok = true
                output_tag="#3"
            end               
            
        elsif mid_loop[0..1]== "CG" and mid_loop[-2..-1]=="AC" and special_A+special_C="AC" and a7+a24="AA" and mid_loop.size  < 7
            #OSTA
            fold_ok = true
            output_tag="#4"
        else  
            if mid_loop.size - check_comp[0]*2 <=1 and special_A+special_C = "AG" and a7+a24 ="AA"
                #Samonella
                fold_ok = true
                output_tag="#5"
            end
        end
    else   
        # puts max_self_complementary(mid_loop)

        if   special_A+special_C=="AG" and a7+a24=="AA"
            if mid_loop.size - check_comp[0]*2 == 0 
                #S9 Pep
               fold_ok = true
               output_tag="#6"
            end
        end 
    end

    unless fold_ok
        unless   rna_full.get_pair(stemI_left_1) == stemI_right_1 and  rna_full.get_pair(stemI_left_2) == stemI_right_2
            warning(warning_on,"#{fasta.header}:obvious global alternative folding\n#{rna_full.folding_str}\n#{rna_full.seq}")
            next
        end
    end


    # test mid_loop
    # when the inner loop is less than 13, which is an arbitary number
    # then test if there is self complementary more than 2 bases, if yes, good to go
    # if no, if the structure has self complementary in the inner side, 
    # if yes, if the structure implicates with ends. no? test less stringent pairs



    unless fold_ok
        if options.has_key?(:s) and options.has_key?(:s)==true
            warning(warning_on,"#{fasta.header}:rejected not special case")
            next
        end
 
        case mid_loop.size
        when 1        
            next
        when 2..7 
            check_comp = max_self_complementary(mid_loop)
            if check_comp[0]>0 and mid_loop.size - check_comp[0]*2 < 2 and check_comp[2]<2                              
                stemII_folding_ok = true
                output_tag<<"_#a"
                #add code here to identify Tr9Tr1 condition
            elsif rna_full.unpaired_between?(special_A_idx,special_C_idx) and rna_full.unpaired_between?(a7_idx,a7_idx+3) 
                #  puts ">#{a7},#{special_A},#{special_C},#{a24},#{mid_loop}"
                if a7=="A" and special_A=="A" and special_C=="A"  and a24=="C" and mid_loop[0..1]=="AG" and mid_loop[-2..-1]=="TG"
                    #case of INVX-19 
                    stemII_folding_ok = true
                    output_tag<<"_#b_1"
                elsif a7=="A" and special_A=="A" and a24=="C" and special_C=="G" and mid_loop[0..1]=="CA" and mid_loop[-2..-1]=="GT"
                    #case of INVX-2
                    stemII_folding_ok = true
                    output_tag<<"_#b_2"
                elsif a7=="A" and a24=="A" and special_A=="C" and special_C=="C" and mid_loop=="AGG"
                    #case INVX-13
                    stemII_folding_ok =true
                    output_tag<<"_#b_4"
                elsif  a7=="A" and a24=="C" and special_A=="A" and special_C=="G" and mid_loop=="CCT"
                    #case INVX-14
                    stemII_folding_ok =true
                    output_tag<<"_#b_5"
                elsif a7=="A" and a24=="C" and special_A=="A" and special_C=="G" and mid_loop[0..1]=="GA" and mid_loop[-2..-1]=="GT"
                    #case INVX_17
                    stemII_folding_ok =true
                    output_tag<<"_#b_6"
                end
            elsif a7=="A" and special_A=="T" and a24=="A" and special_C=="G" and is_paired?(fasta.seq[special_A_idx+1], fasta.seq[special_C_idx-1]) and rna_full.unpaired_between?(special_A_idx+1, special_C_idx) 
                output_tag<<"_#b_7"
                stemII_folding_ok = true
            
            elsif check_comp[0]>1
                output_tag<<"_#b_8"
                stemII_folding_ok = true
            else
                warning(warning_on,"#{fasta.header}:rejected short mid_loop pairing 2..7")
                next                
            end

            # test if there is some complementary with other part of the RNA
        when 8..16

            if  is_paired?(mid_loop[0],mid_loop[-1]) and mid_loop.size - check_comp[0]*2 < mid_loop.size/2 and check_comp[2]<2
                #Gordi case NZ_RKMH01000011.1
                stemII_folding_ok = true
                output_tag<<"_#c"
            
            elsif is_paired?(mid_loop[0],mid_loop[-1]) and rna_full.paired?(special_A_idx+1,special_C_idx-1)        
                stemII_folding_ok = true
                output_tag<<"_#d"
            elsif rna_full.unpaired_between?(special_A_idx,special_C_idx) and rna_full.unpaired_between?(a7_idx,a7_idx+3)
                if  a7=="A" and a24=="A" and special_A=="A" and special_C=="G" and mid_loop[0..3]=="GAAT" and mid_loop[-4..-1]=="AACA"
                    #case INVX_16
                    stemII_folding_ok = true
                    output_tag<<"_#e_1"
                elsif a7=="A" and special_A=="C" and a24=="A" and special_C=="T" and mid_loop[0..1]=="AC" and mid_loop[-2..-1]=="GA"
                    # case of INV-X6
                    stemII_folding_ok = true
                    output_tag<<"_#e_2" 
                end
            else
                unless is_paired?(mid_loop[0],mid_loop[-1])
                    warning(warning_on,"#{fasta.header}:rejected no stemII first base pairing ")
                    next
                end

                count_in_mid_fold = rna_mid.pair_count(0,mid_loop.size-1)
                count_in_full_fold = rna_full.pair_count(special_A_idx+1,special_C_idx-1)
                # puts count_in_mid_fold, count_in_full_fold,rna_mid.left_unfolded,rna_mid.right_unfolded,rna_mid.left_unfolded,rna_mid.right_unfolded
                if count_in_full_fold>=3 and count_in_full_fold-count_in_mid_fold <2 
                    if (rna_mid.left_unfolded-rna_mid.right_unfolded).abs <2 and [rna_mid.left_unfolded,rna_mid.right_unfolded].max <4
                        stemII_folding_ok = true  
                        output_tag<<"_#f"
                        # long S9
                    end                                 
                end
           
            end


            
        when 17..50
            # folding information required    
            unless is_paired?(mid_loop[0],mid_loop[-1])
                warning(warning_on,"#{fasta.header}:rejected no stemII ")
                next
            end
    
            if rna_full.paired?(special_A_idx+1,special_C_idx-1)        
                stemII_folding_ok = true
                output_tag<<"_#g"
               
            else
               
                count_in_mid_fold = rna_mid.pair_count(0,mid_loop.size-1)
                count_in_full_fold = rna_full.pair_count(special_A_idx+1,special_C_idx-1)
                if count_in_full_fold>7 and count_in_full_fold-count_in_mid_fold <3 
                    if (rna_mid.left_unfolded-rna_mid.right_unfolded).abs <2 and [rna_mid.left_unfolded,rna_mid.right_unfolded].max <4
                        if ["CG","GC"].include?(fasta.seq[special_A_idx+1]+fasta.seq[special_C_idx-1])
                            stemII_folding_ok = true 
                            output_tag<<"_#h"
                             # invx12 case
                        end
                    end                                 
                end
           
            end
            # other tests

        else            
            warning(warning_on,"#{fasta.header}:rejected mid_loop length more than 50 or < 2")
            next
        end
        # end of mid_loop test


        if is_paired?("CG"+special_A,(mid_loop[-2..-1]+special_C).reverse)
            # case of BAC_Methylobacterium
            # CGA ... TCG
            warning(warning_on,"rejection:#{fasta.header},case XM_015126877")
            stemII_folding_ok=false
        end

        
        if reverse_complement(special_A+mid_loop[0]) == special_C+a24
            # depletion of free A
            warning(warning_on,"rejection:#{fasta.header},case spcecial_A and mid_loop0:spcial_C+a24 paring")
            stemII_folding_ok=false
        end
        
       
        unless stemII_folding_ok
            warning(warning_on,"#{fasta.header}:rejected no stemII not OK ")
            next
        end

        output_tag="#7"+output_tag
     
        
    end

    drawing<<"."*(seq_head.size-stemI_len+1) # 
    drawing<<"("*stemI_len
    if a7=='A'
        drawing<<"*"
    else
        drawing<<"&"
    end
    drawing<< "*****"
    if special_A=="A"
        drawing<<"~"
    else
        drawing<<"-"
    end

    inner_stem_len= rev_comp_length(mid_loop,mid_loop,false)
    if inner_stem_len*2 >mid_loop.size
        drawing<<"."*mid_loop.size
    else
        drawing<<'('*inner_stem_len+'.'*(mid_loop.size-2*inner_stem_len)+')'*inner_stem_len
    end
    if special_C=="C"
        drawing<<"~"
    else
        drawing<<"-"
    end

    drawing<<"*****"
    drawing<<")"*stemI_len
    drawing<<"."*(seq_tail.size-stemI_len+1)# 

    if options.has_key?(:r) and options[:r]==true
        next
    end

    puts fasta.header
    if options[:i]==true
        puts output_tag
    end
    puts fasta.seq
    if options.has_key?(:d) and options[:d]
        puts drawing
        # if warning_on
        #     puts "tripleII_A: #{a24}"
        #     puts "Special_A:#{special_A}"
        #     puts "Special_C:#{special_C}"
        #     puts "mid_loop:#{mid_loop}"
        # end  
    end

end

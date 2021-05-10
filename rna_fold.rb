# Author: Hengyi Jiang <hengyi.jiang@gmail.com>
# Purpose: library for wrapping Vienna RNAfold program and fold result parsing
# Requirement: Vienna RNA Package

require "fileutils"
require "digest"
require_relative("fasta.rb")

class Fold_parser
	# this include the 
	def initialize
		@folding_str=""
		@pairs=[]
	end
	def continous(x,y)
	 	return ( x==y-1 or x== y+1 )
	end 

	def parse(str)
		# to generate the @pair which is an array of array each element representing
		# [pair_left, pair_right, pair_lenth]
		#      |            |                       
		#  ....(((((....)))))....     5  

		@pairs.clear
		@folding_str=str
		front_stack=[]
		back_stack=[]
		@pair_stack =[]
		# ..(((..)))   pair_stack=[[4,7],[3,8],[2,9]]
		layer=0
		last_tag=""
		for i in 0..(str.size-1)
			tag = str[i]
			if tag=="("
				front_stack << i
			elsif tag==")"
				@pair_stack <<[front_stack.pop,i]
			end	
		end
		@pair_stack.sort! {|x,y| x[0]<=>y[0]}
		#p @pair_stack
		last_left=-1
		last_right=-1
		@pair_stack.each do |pp|
			if continous(last_left, pp[0])  and continous(last_right,pp[1])
				 # do nothing
			else
				if last_left ==-1
					@pairs<<[pp[0],pp[1],1]
				elsif
					@pairs[-1][1]=last_right
					@pairs[-1][2]=last_left-@pairs[-1][0]+1
					@pairs<<[pp[0],pp[1],1]
				end
			end

			last_left=pp[0]
			last_right=pp[1]
		end
		if @pair_stack.size>1
			@pairs[-1][1]=last_right
			@pairs[-1][2]=last_left-@pairs[-1][0]+1
		end
		# p @pairs
	end
	attr_accessor :pairs, :pair_stack
end 

class Fold_file_parser
	def initialize(file_name)
		begin
			@fh = File.open(file_name)
		rescue
			puts ("file open failed #{file_name}")
		end 
	end

	def close
		@fh.close
	end
	def next
		n=1
		while n<4 and line=@fh.gets 
			t = n%3
			case t
			when 1
				header=line.chomp
			when 2
				seq=line.chomp
				
			when 0
				if m=line.match(/([\.\(\)]+)\s+\(\s?(.*)\)$/)
					folding_str=m[1]
					folding_energy=m[2].to_f
					rna=Folded_RNA.new(seq)
					rna.header=header
					rna.folding_str=folding_str
					
					fp=Fold_parser.new
					fp.parse(folding_str)
					rna.pairs=fp.pairs
					rna.pair_stack=fp.pair_stack
					rna.folding_energy=folding_energy
				else
					
					raise("File format problem")
				end

			end
			n +=1
		end
		rna 	
	end
end

class Folded_RNA
	def initialize(seq)
		@seq=seq
		@folding_str=nil
		@folding_energy=nil
		@pairs=[]
		@pair_stack=[]
		@header=""
		@is_closed
	end	
	def fold!(temperature=37)
		file_prefix="rnafold_#{Time.new.to_i}"
		f=File.new("#{file_prefix}.tmp","w")
		f.puts @seq
		f.close
		cmd ="RNAfold  --noPS -i #{file_prefix}.tmp -T #{temperature.to_s} > #{file_prefix}.folding"
		ret=`#{cmd}`
		unless ret.empty?
			raise "RNAfold error[#{cmd}]: #{ret}"
		end
		if File.exist?( "#{file_prefix}.tmp")
			FileUtils.rm("#{file_prefix}.tmp")
		end

		f=File.open("#{file_prefix}.folding")
		f.readline()
		folding=f.readline
		f.close

		if m=folding.match(/([\.\(\)]+)\s+\(\s?(.*)\)$/)
			@folding_str=m[1]
			@folding_energy=m[2].to_f
		end

		if @folding_str !=nil
			parser=Fold_parser.new
			parser.parse(@folding_str)
			@pairs=parser.pairs
			@pair_stack=parser.pair_stack
		end 
		
		#puts folding
		if File.exist?("#{file_prefix}.folding")
			FileUtils.rm("#{file_prefix}.folding")
		end
	end

	def left_unfolded
		if @folding_str==nil
			fold
		end
		if m=@folding_str.match(/^(\.+)/)
			m[1].size
		else
			0
		end
	end
	def right_unfolded
		if @folding_str==nil
			fold
		end
		if m=@folding_str.match(/(\.+)$/)
			m[1].size
		else
			0
		end
	end

	def is_stem_loop?
		if @folding_str==nil
			fold
		end
		if @pairs.size ==1
			true
		else
			false
		end
	end
	def is_closed?
	#  the definition of closed structure means the ends of the RNA pairs with each other
    #  .....((((..........((((...))))......))))...

		if @is_closed==nil
			if @folding_str==nil
				fold
			end
			if @pairs.size ==1
				@is_closed=true
				return true
			elsif @pairs.size>1
				r=@pairs[0][1]
				for i in (1..(@pairs.size-1))
					if r < @pairs[i][1]
						not_closed=true
						break
					end
				end 
				if not_closed
					@is_closed=false
					return false
				else
					@is_closed=true
					return true
				end
			else
				return false
			end
		else
			return @is_closed 
		end
	end

	def closed_stem_len
		if is_closed?
			return @pairs[0][2]
		else
			return 0
		end
	end


	def not_structured?
		if @folding_str==nil
			fold
		end
		if @folding_str.match(/^\.*$/)
			return true
		else
			return false
		end
	end

	def stem_len
		if is_stem_loop?
			@pairs[0][2]
		else
			0
		end
	end

	def loop_len
		if is_stem_loop?
			@pairs[0][1]-@pairs[0][0]-@pairs[0][2]
		else
			0
		end
	end
	def pairing_tag(left,right)
		left_=left-1
		right_=right-1
		v =""
		r=[]
		for i in left_.upto(right_)
			for pp in @pair_stack 
				if i==pp[0]
					r<<pp[1]
				elsif i==pp[1]
					r<<pp[0]
				else
					r<< -1
				end
			end
		end
		for i in (0..(@seq.size-1) )
			if (left_..right_).include?(i)
				v << "-"
			elsif r.include?(i)
				v << "*"
			else 
				v << " "
			end
		end
		return v
	end

	def folding_str
		if @folding_str==nil
			fold
		end
		@folding_str
	end

	def get_pair(n)
		if @pair_stack ==nil
			return nil
		end
		@pair_stack.each do |pr|
			if pr[0]==n
				return pr[1]
			elsif pr[1]==n
				return pr[0]
			end			
		end
		return -1
	end

	def get_stretch(n)
		if @pair_stack ==nil
			return []
		end
		@pairs.each do |pr|
			if (pr[0]..(pr[0]+pr[2]-1)).include?(n)
				return [pr[0],(pr[0]+pr[2]-1)]
			elsif ((pr[1]-pr[2]+1)..pr[1]).include?(n)
				return [pr[1]-pr[2]+1,pr[1]]
			end
		end
		return []
	end

	def unpaired?(index_array)
		if index_array.is_a?(Array)
			index_array.each do |i|
				if get_pair(i) != -1
					return false
				end
			end
			return true
		end
		return get_pair(index_array)==-1
		
	end

	def unpaired_between?(a,b)
		if a>b
			a,b = b,a
		end
		l=a.upto(b).map{|n| n}
		return unpaired?(l)
	end
	

	def paired?(a,b)
		a_pair = get_pair(a)
		if a_pair == nil
			return false
		end
		if a_pair==b
			return true
		end
	end

	def pair_count(a,b)
		ar=[]
		@pair_stack.each do |s|
			if (a..b).include?(s[0]) or (a..b).include?(s[1])
				if !ar.include?(s[0].to_s+s[1].to_s)
					ar << s[0].to_s+s[1].to_s
				end
			end
		end
		return ar.size
	end


	def next_pair(idx)

		tmp=[]
		#case 1

		if @pair_stack.empty?
			return nil
		end

		#.....(((....(((...))..)))).X..
		#                           |--idx
		if @pair_stack[-1][1] <= idx
			return nil
		end

		@pair_stack.each do |pr|
			tmp = pr
			if pr[0] <= idx
				next
			else					
				break
			end
		end

		# case 2
		#.....(((.X..(((...))..))))
		#         |--idx
		if tmp[0]>idx
			return tmp
		end
		
		#case 2
		#.....(((...(((..)).x.))))
		#                   |--idx
		@pair_stack.reverse_each do |pr|
			tmp = pr
			if pr[1] <= idx
				next
			else					
				break
			end
		end
		return tmp
	end


	attr_accessor :seq, :header,:folding_str,:folding_energy,:pairs, :pair_stack
end


class Folded_collection
	include Enumerable

	def initialize
		@folded_RNAs=[]
		@index_hash={}
	end
	def fold_fasta_collection!(ft_collection,temperature=37)
		unless ft_collection.is_a?(Fasta_collection)
			raise "Folded_collection::fold_fasta_collection:Fasta_collection is expected"
		end

		file_prefix="rnafold_#{Time.new.to_i}"
		tmp_file=file_prefix+"_tmp.fasta"
		begin
			fh=File.open(tmp_file,"w")
		rescue
			#retry			
			file_prefix="rnafold_#{Time.new.to_i}"
			if !File.exist?(tmp_file)
				fh=File.open(tmp_file,"w")
			else
				raise "can't solve the problem of repeated file name:#{tmp_file}"
			end
		end

		ft_collection.each do |fasta|
			fh.puts fasta.header
			fh.puts fasta.seq			
		end
		fh.close

		fold_cmd ="RNAfold  --noPS -T #{temperature} -i #{tmp_file} > #{file_prefix}.folding"
		ret=`#{fold_cmd}`
		unless ret.empty?
			raise "RNAfold error[#{cmd}]: #{ret}"
		end

		unless File.exist?("#{file_prefix}.folding")
			raise "file #{file_prefix}.folding not exists"		
		end
		parser=Fold_file_parser.new("#{file_prefix}.folding")
		i=0
		while rna=parser.next
			@folded_RNAs<<rna
			md5 = Digest::MD5.hexdigest(rna.header)
			@index_hash[md5]=i
			i +=1
		end
		if File.exist?( tmp_file)
			FileUtils.rm(tmp_file)
		end

		if File.exist?( "#{file_prefix}.folding")
			FileUtils.rm("#{file_prefix}.folding")
		end

		true
	end

	def each
		@folded_RNAs.each do |obj|
			yield obj
		end
	end

	def get(header)
		md5 = Digest::MD5.hexdigest(header)
		if @index_hash.has_key?(md5)
			return @folded_RNAs[@index_hash[md5]]
		else
			nil
		end
	end
	def size
		@folded_RNAs.size
	end

end

# #testing codes
# puts "testing Fold_parser"
# test_string="...(((((..)))...(((..))).))..."
# parser=Fold_parser.new
# parser.parse(test_string)
# test_string="(((((..............((((......))))(((((((((((...)))))))))..))(((((((((.............))))))))))))))"
# parser.parse(test_string)
# p parser

# puts "testing Folded_RNA"
# rna= Folded_RNA.new("ACCTATACGGTATCCGCATTCCAATGC")
# rna.fold
# p rna
# puts rna.left_unfolded
# puts rna.right_unfolded
# if  rna.is_stem_loop?
# 	puts "stem_len:#{rna.stem_len}"
# 	puts "loop_len:#{rna.loop_len}"
# end
# puts rna.seq
# puts rna.folding_str
# puts rna.pairing_tag(3,5)

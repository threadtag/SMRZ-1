# Author: Hengyi Jiang <hengyi.jiang@gmail.com>
# liberary of functions commmonly used to handle sequences

require_relative 'fasta.rb'
def reverse_complement(seq)
	replace_map ={'A'=>"T","C"=>"G","T"=>"A","G"=>"C","a"=>'t','c'=>'g','t'=>'a','g'=>'c',"\n"=>""}
	r=seq.reverse
 	r.gsub(Regexp.union(replace_map.keys), replace_map)
end


def count_base(in_seq)
# input sequence string, RNA or DNA
# return hash of Base count
# N any Char, X other Non-DNA or RNA char

	r={"A"=>0,"T"=>0,"C"=>0,"G"=>0,"U"=>0,"X"=>0,"N"=>0}
	in_seq.each_char do |s|
		if s.match(/[ATCGU]/i)
			r[s.upcase]=r[s.upcase]+1
		else
			r['X']=r['X']+1
		end
		r['N'] = r['N']+1
	end
	return r
end

def self_comp_check(in_seq)
# input : normal sequence string, RNA or DNA
# return: array 
# [match_i[max_t],match_span[max_t],match_j[max_t],match_seq_l[max_t],match_seq_r[max_t],match_seq_score[max_t]]
# or [0,0,0]
#    _span__                    _span__
#	 i   span-1             j-span+1  j
#    |     |                    |     |
#	[a][b][c][d][e][f][g][h][i][j][k][l]
#   i scan start from left, with min span of 2, 
#   j scan start from right, with min span equal to i-span

	len=in_seq.size
	match_flag =0
	match_i=[]
	match_span=[]
	match_j=[]
	match_seq_l=[]
	match_seq_r=[]
	match_seq_score=[]
	for i in 0.upto(len-3)
		for span in 2.upto(len/2)
			for j in (len-1).downto(i+2*(span-1)+1)
				tester_l=in_seq[i..(i+span-1)]
				tester_r=in_seq[(j-span+1)..j]
				
				if is_paired?(tester_l.upcase, tester_r.upcase.reverse, true) # allow wobble G-T pair
					match_seq_l<<tester_l.upcase
					match_seq_r<<tester_r.upcase
					tester_composition=count_base(tester_l.upcase)
					match_seq_score<<tester_composition['A']*2+tester_composition['T']*2+tester_composition['C']*3+tester_composition['G']*3
					match_flag =1
					match_i<<i
					match_span<<span
					match_j<<j
				end	
			end

			if match_flag ==1
				match_flag=0				
				next; #next span for match,			
			else 
				break;# no further test for other span, goto next i 
			end
		end
	end

	max_score=0
	max_t=0
	max_span=0
	# find the max score match, if score is equal, choose the longer one
	for t in 0.upto(match_i.size-1)
		if match_seq_score[t]>max_score
			max_score = match_seq_score[t]
			max_t =t
			max_span=match_span[t]
		elsif match_seq_score[t]==max_score
			if match_span[t] > max_span
				max_t =t
				max_span=match_span[t]
			end			
		end
	end
	if !match_i.empty?
		#return [match_i[max_t],match_span[max_t],match_j[max_t],match_seq_l[max_t],match_seq_r[max_t],match_seq_score[max_t]]
		return {:ok=>true,:left_i=>match_i[max_t],:span=>match_span[max_t],:right_i=>match_j[max_t]-match_span[max_t]+1,
				:left_seq=>match_seq_l[max_t],:righ_seq=>match_seq_r[max_t],:score=>match_seq_score[max_t],
				:unpaired_left=>match_i[max_t],:unpaired_right =>len-match_j[max_t]-1,
				:loop_len=>match_j[max_t]-match_i[max_t]-2*match_span[max_t]+1  #j-i-2*span+1
			   }		
				#    _span__                    _span__
				#	 i   span-1             j-span+1  j
				#    |     |                    |     |
				#	[a][b][c][d][e][f][g][h][i][j][k][l]
	else
		return {:ok=>false}
	end
end


def max_self_complementary(seq)
	half_len = seq.size/2-1
	max_continious = 0 # disrupting continuious
	last_continious = 0
	m=0
	n=0
	t=0
	0.upto(half_len).each do |l|
		if is_complement?(seq[l],seq[-(l+1)],true)
			m +=1
			t=0
			
			if max_continious<last_continious
				max_continious=last_continious
			end
			last_continious=0
		else

			n +=1
			if t ==1
				last_continious +=1
			else
				if n ==1
					last_continious=1
				end			
			end
			t=1			
		end
	end

	return [m, n, max_continious]
end 
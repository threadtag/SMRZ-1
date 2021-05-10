# Author: Hengyi Jiang <hengyi.jiang@gmail.com>
class Header_parser

    def initialize(h)
        @header=h
        parse
    end

    def parse
        if @m=@header.match(/>(\S+)\s(.+)_(\w\w)_(\d+)_(\d+)/)
            @ok=true
            return true
        else
            @ok=false
            return false
        end
    end

    def id
        if @ok
            return @m[1]
        else
            return @header
        end
    end

    def genome
        if @ok
            return @m[2]
        else
            return @header
        end
    end

    def genome_short
        if @ok
            if @m[2].match(/strain/i)
                i= @m[2].index("strain")
                return @m[2][0..i-2]
            else
                return @m[2].split(/\s+/)[0..1].join(" ")
            end
        else
            return @header
        end
    end

    def start
        if @ok
            return @m[4].to_i
        else
            return -1
        end
    end

    def stop
        if @m
            return @m[5].to_i
        else
            return -1
        end
    end

    def plus_strand?
        if @m
            if @m[3]=='fw'
                return true
            else
                return false
            end

        else
            return false
        end
    end

    def minuse_strand?
        if @m
            if @m[3]=='rv'
                return true
            else
                return false
            end
        else
            return false
        end
    end
    def which_strand
        if plus_strand?
            return "+"
        elsif minuse_strand?
            return "-"
        else
            "?"
        end
    end

    def ok?
        @ok
    end
end
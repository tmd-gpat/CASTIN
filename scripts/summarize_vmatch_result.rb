refseq_table = []
refseq_gene_table = {}
entries = {}

if ARGV.length < 4
  exit
end
collect_d = ARGV[3].include? 'd'
collect_p = ARGV[3].include? 'p'

STDERR.puts "loading reference fasta."
open(ARGV[0]) { |file|
  while l = file.gets do
    if l[0] == '>' then
      refseq = l.split()[0].slice(1,100)
      refseq_table << refseq
      entries[refseq] = Array.new
    end
  end
}

STDERR.puts "loading refLink."
open(ARGV[1]) { |file|
  while l = file.gets do
    row = l.split("\t")
    refseq = row[2]
    gene_id = row[6].to_i
    refseq_gene_table[refseq] = gene_id
  end
}

STDERR.puts "loading vmatch output."
open(ARGV[2]) { |file|
  n = 0
  while l = file.gets do
    next if l[0] == '#'

    n+=1
    STDERR.puts "processed #{n} lines." if n % 1000000 == 0

    row = l.split
    refseq1 = refseq_table[row[1].to_i]
    refseq2 = refseq_table[row[5].to_i]
    length1 = row[0].to_i
    length2 = row[4].to_i
    position1 = row[2].to_i
    position2 = row[6].to_i
    direction = row[3]
    next if direction == 'D' && !collect_d
    next if direction == 'P' && !collect_p

    next if refseq_gene_table[refseq1].nil? || refseq_gene_table[refseq2].nil?

    # DEBUG: output refseq-refseq
    # puts "#{refseq1}\t#{refseq2}\t#{position1}\t#{length1}\t#{position2}\t#{length2}"

    # discard splice variant
    if refseq_gene_table[refseq1] != refseq_gene_table[refseq2]
      new_entry = true
      i = 0
      while i < entries[refseq1].length do
        entry = entries[refseq1][i]
        if position1 <= entry[:from] && entry[:to] <= position1 + length1 - 1
          entries[refseq1].delete_at i
          i -= 1
        elsif entry[:from] <= position1 && position1 + length1 - 1 <= entry[:to]
          new_entry = false
        end
        i += 1
      end
      if new_entry
        entries[refseq1] << {
          refseq: refseq2,
          from: position1,
          to: position1 + length1 - 1,
          from_target: position2,
          to_target: position2 + length2 - 1,
        }
      end

      new_entry = true
      i = 0
      while i < entries[refseq2].length do
        entry = entries[refseq2][i]
        if position2 <= entry[:from] && entry[:to] <= position2 + length2 - 1
          entries[refseq2].delete_at i
          i -= 1
        elsif entry[:from] <= position2 && position2 + length2 - 1 <= entry[:to]
          new_entry = false
        end
        i += 1
      end
      if new_entry
        entries[refseq2] << {
          refseq: refseq1,
          from: position2,
          to: position2 + length2 - 1,
          from_target: position1,
          to_target: position1 + length1 - 1,
        }
      end
    end
  end
}

refseq_table.each { |refseq|
  if entries[refseq].size > 0
    puts (">" + refseq)
    entries[refseq].sort {|a,b|
      a[:from] <=> b[:from]
    }.each {|entry|
      puts "#{entry[:refseq]}\t#{entry[:from]}\t#{entry[:to]}\t#{entry[:from_target]}\t#{entry[:to_target]}"
    }
  end
}


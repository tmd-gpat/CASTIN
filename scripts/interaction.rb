# Symbol_human, Symbol_mouseを入力として, 各種KEGGHPRD_*ファイルを出力する
#
# $ ruby scripts/interaction.rb ../output/ExpID-XYZ/ parameters/

require 'set'

## constants
deprecated_entrez_ids = [
#    6081,   # SNORA73B (human)
#    404338, # Olfr1434 (mouse)
]
interaction_path = "KEGG/curated_HPRD_KEGG_20161110.txt"

## start of class definition
class Gene
    attr_accessor :symbol, :taxonomy, :entrez_id    # refLink
    attr_accessor :homologue    # homologene
    attr_accessor :expression, :raw_count

    def initialize(symbol, entrez_id, taxonomy)
        @symbol = symbol
        @entrez_id = entrez_id
        @taxonomy = taxonomy
    end
end
class Interaction
    attr_accessor :ligand_symbol, :receptor_symbol, :type, :id, :kegg
    attr_accessor :ligand_human, :ligand_mouse, :receptor_human, :receptor_mouse

    def initialize(ligand_symbol, receptor_symbol, type, id, kegg,
                   ligand_human, ligand_mouse, receptor_human, receptor_mouse)
        @ligand_symbol = ligand_symbol
        @receptor_symbol = receptor_symbol
        @type = type
        @id = id
        @kegg = kegg
        @ligand_human = ligand_human
        @ligand_mouse = ligand_mouse
        @receptor_human = receptor_human
        @receptor_mouse = receptor_mouse
    end

    def url
        "http://www.genome.jp/kegg/pathway/hsa/hsa" + @id + ".html";
    end
end
## end of class definition

## start of main process
STDERR.puts("procesing directory: #{ARGV[0]}")

## load parameter files
# load refseq_id list for each taxonomy
refseq_ids_human = Set.new
refseq_ids_mouse = Set.new

open(ARGV[1] + "/hg38/refNames.txt"){ |file|
    while l = file.gets do
        refseq_ids_human.add(l.chomp)
    end
}
open(ARGV[1] + "/mm10/refNames.txt"){ |file|
    while l = file.gets do
        refseq_ids_mouse.add(l.chomp)
    end
}
STDERR.puts("#{refseq_ids_human.length} human refseqs & #{refseq_ids_mouse.length} mouse refseqs from refNames files are loaded.")

# load refLink to relate symbol, entrez_id and taxonomy
genes_human = Hash.new
genes_mouse = Hash.new
genes_human_by_entrez = Hash.new
genes_mouse_by_entrez = Hash.new

open(ARGV[1] + "/refLink.txt"){ |file|
    while l = file.gets do
        row = l.chomp.split("\t")
        symbol = row[0]
        refseq_id = row[2]
        entrez_id = row[6].to_i
        next if (deprecated_entrez_ids.include? entrez_id)
        if refseq_ids_human.include? refseq_id then
            if genes_human.include? symbol and genes_human[symbol].entrez_id != entrez_id then
                STDERR.puts("duplicated gene: #{symbol} : #{entrez_id} & #{genes_human[symbol].entrez_id} (human)")
            else
                gene = Gene.new(symbol, entrez_id, "human")
                genes_human[symbol] = gene
                genes_human_by_entrez[entrez_id] = gene
            end
        elsif refseq_ids_mouse.include? refseq_id then
            if genes_mouse.include? symbol and genes_mouse[symbol].entrez_id != entrez_id then
                STDERR.puts("duplicated gene: #{symbol} : #{entrez_id} & #{genes_mouse[symbol].entrez_id} (mouse)")
            else
                gene = Gene.new(symbol, entrez_id, "mouse")
                genes_mouse[symbol] = gene
                genes_mouse_by_entrez[entrez_id] = gene
            end
        end
    end
    STDERR.puts("#{genes_human.size} human genes & #{genes_mouse.size} mouse genes from refLink file are loaded.")
}

# load homologene.data to relate human/mouse symbols
open(ARGV[1] + "/homologene.data"){ |file|
    current_homologene_id = -1
    human_entrez_id = -1
    mouse_entrez_id = -1
    homologue_count = 0
    while l = file.gets do
        row = l.chomp.split("\t")
        homologene_id = row[0].to_i
        tax_id = row[1].to_i
        entrez_id = row[2].to_i
        symbol = row[3]

        if current_homologene_id != homologene_id then
            human_entrez_id = -1
            mouse_entrez_id = -1
            current_homologene_id = homologene_id
        end
        if tax_id == 9606 then
            human_entrez_id = entrez_id
            if mouse_entrez_id != -1 then
                if genes_human_by_entrez.include? human_entrez_id and
                   genes_mouse_by_entrez.include? mouse_entrez_id then
                    genes_human_by_entrez[human_entrez_id].homologue = mouse_entrez_id
                    genes_mouse_by_entrez[mouse_entrez_id].homologue = human_entrez_id
                    homologue_count += 1
                end
            end
        elsif tax_id == 10090 then
            mouse_entrez_id = entrez_id
            if human_entrez_id != -1 then
                if genes_human_by_entrez.include? human_entrez_id and
                   genes_mouse_by_entrez.include? mouse_entrez_id then
                    genes_human_by_entrez[human_entrez_id].homologue = mouse_entrez_id
                    genes_mouse_by_entrez[mouse_entrez_id].homologue = human_entrez_id
                    homologue_count += 1
                end
            end
        end
    end
    STDERR.puts("#{homologue_count} homologue entries are loaded.")
}

# load interactions
interactions = []
open(ARGV[1] + interaction_path){ |file|
    file.gets
    while l = file.gets do
        row = l.scrub.chomp.split("\t")
        interaction_id = row[0].to_i
        ligand = row[1]
        receptor = row[2]
        kegg = "NA"
        id = "NA"
        if (row[10].split("_").size >= 2) then
            id = row[10].split("_")[0]
            kegg = row[10].split("_")[1]
        end
        type = row[11]

        # search human/mouse entrez_ids
        ligand_human = nil
        ligand_mouse = nil
        receptor_human = nil
        receptor_mouse = nil

        if genes_human[ligand] != nil then
            ligand_human = genes_human[ligand]
            ligand_mouse = genes_mouse_by_entrez[genes_human[ligand].homologue]
        end
        if genes_human[receptor] != nil then
            receptor_human = genes_human[receptor]
            receptor_mouse = genes_mouse_by_entrez[genes_human[receptor].homologue]
        end

        interaction = Interaction.new(
            ligand, receptor, type, id, kegg,
            ligand_human, ligand_mouse, receptor_human, receptor_mouse)
        interactions << interaction
    end
    STDERR.puts("#{interactions.size} interactions loaded.")
}

## load Symbol files
open(ARGV[0] + "/Symbol_cancer.txt"){ |file|
    file.gets
    while l = file.gets do
        row = l.chomp.split("\t")
        symbol = row[0]
        refseq = row[1]
        raw = row[2].to_i
        expression = row[4].to_f

        if genes_human[symbol] != nil then
            genes_human[symbol].raw_count = raw
            genes_human[symbol].expression = expression
        else
            STDERR.puts("unknown gene : #{symbol} (human)")
        end
    end
}
open(ARGV[0] + "/Symbol_stroma.txt"){ |file|
    file.gets
    while l = file.gets do
        row = l.chomp.split("\t")
        symbol = row[0]
        refseq = row[1]
        raw = row[2].to_i
        expression = row[4].to_f

        if genes_mouse[symbol] != nil then
            genes_mouse[symbol].raw_count = raw
            genes_mouse[symbol].expression = expression
        else
            STDERR.puts("unknown gene : #{symbol} (mouse)")
        end
    end
}
STDERR.puts("loaded symbol files.")

## analysis

# KEGGHPRD_result
open(ARGV[0] + "/KEGGHPRD_result_reanalysis.txt", "w"){ |file|
    # header line
    file.puts("ligand\treceptor\tcount(ligand)\tcount(receptor)\traw count(ligand)\traw count(receptor)\taverage\tligand ratio\tcount of the other ligands\treceptor ratio\tinteraction type\tpathway\tlink")
    interactions.each{ |inter|
        # count sum of other ligand expression
        sum_of_other_ligand_expression = 0.0
        interactions.each { |inter2|
            # symbols
            if inter.receptor_symbol == inter2.receptor_symbol then
                if inter2.ligand_human != nil then
                    sum_of_other_ligand_expression += inter2.ligand_human.expression
                end
                if inter2.ligand_mouse != nil then
                    sum_of_other_ligand_expression += inter2.ligand_mouse.expression
                end
            end
        }

        # expressions
        file.print("#{inter.ligand_symbol}(cancer#{inter.ligand_human == nil ? "/NA" : ""})\t")
        file.print("#{inter.receptor_symbol}(stroma#{inter.receptor_mouse == nil ? "/NA" : ""})\t")
        file.print("#{inter.ligand_human != nil ? inter.ligand_human.expression : "NA"}\t")
        file.print("#{inter.receptor_mouse != nil ? inter.receptor_mouse.expression : "NA"}\t")
        file.print("#{inter.ligand_human != nil ? inter.ligand_human.raw_count : "NA"}\t")
        file.print("#{inter.receptor_mouse != nil ? inter.receptor_mouse.raw_count : "NA"}\t")

        # ligand (human) / receptor (mouse)
        if inter.ligand_human != nil && inter.receptor_mouse != nil then
            # average
            average = Math.sqrt(inter.ligand_human.expression * inter.receptor_mouse.expression)
            file.print("#{average}\t");
            # lig human/all ratio
            lig_expression_all = 0
            lig_expression_all += inter.ligand_human.expression if inter.ligand_human != nil
            lig_expression_all += inter.ligand_mouse.expression if inter.ligand_mouse != nil
            if lig_expression_all > 0 then
                lig_ratio_human = inter.ligand_human.expression / lig_expression_all
                file.print("#{lig_ratio_human}\t")
            else
                file.print("NA\t")
            end
            # lig/other lig ratio to the receptor
            if sum_of_other_ligand_expression > 0 then
                file.print("#{lig_expression_all / sum_of_other_ligand_expression}\t")
            else
                file.print("0.0\t")
            end
            # rec human/all ratio
            rec_expression_all = 0
            rec_expression_all += inter.receptor_human.expression if inter.receptor_human != nil
            rec_expression_all += inter.receptor_mouse.expression if inter.receptor_mouse != nil
            if rec_expression_all > 0 then
                rec_ratio_mouse = inter.receptor_mouse.expression / rec_expression_all
                file.print("#{rec_ratio_mouse}\t")
            else
                file.print("NA\t")
            end
        else
            file.print("NA\tNA\tNA\tNA\t")
        end
        file.print("#{inter.type}\t#{inter.kegg}\t#{inter.url}\t\n")

        # count sum of other ligand expression

        # expressions
        file.print("#{inter.ligand_symbol}(stroma#{inter.ligand_mouse == nil ? "/NA" : ""})\t")
        file.print("#{inter.receptor_symbol}(cancer#{inter.receptor_human == nil ? "/NA" : ""})\t")
        file.print("#{inter.ligand_mouse != nil ? inter.ligand_mouse.expression : "NA"}\t")
        file.print("#{inter.receptor_human != nil ? inter.receptor_human.expression : "NA"}\t")
        file.print("#{inter.ligand_mouse != nil ? inter.ligand_mouse.raw_count : "NA"}\t")
        file.print("#{inter.receptor_human != nil ? inter.receptor_human.raw_count : "NA"}\t")

        # ligand (mouse) / receptor (human)
        if inter.ligand_mouse != nil && inter.receptor_human != nil then
            # average
            average = Math.sqrt(inter.ligand_mouse.expression * inter.receptor_human.expression)
            file.print("#{average}\t");
            # lig human/all ratio
            lig_expression_all = 0
            lig_expression_all += inter.ligand_human.expression if inter.ligand_human != nil
            lig_expression_all += inter.ligand_mouse.expression if inter.ligand_mouse != nil
            if lig_expression_all > 0 then
                lig_ratio_mouse = inter.ligand_mouse.expression / lig_expression_all
                file.print("#{lig_ratio_mouse}\t")
            else
                file.print("NA\t")
            end
            # lig/other lig ratio to the receptor
            if sum_of_other_ligand_expression > 0 then
                file.print("#{lig_expression_all / sum_of_other_ligand_expression}\t")
            else
                file.print("0.0\t")
            end
            # rec human/all ratio
            rec_expression_all = 0
            rec_expression_all += inter.receptor_human.expression if inter.receptor_human != nil
            rec_expression_all += inter.receptor_mouse.expression if inter.receptor_mouse != nil
            if rec_expression_all > 0 then
                rec_ratio_human = inter.receptor_human.expression / rec_expression_all
                file.print("#{rec_ratio_human}\t")
            else
                file.print("NA\t")
            end
        else
            file.print("NA\tNA\tNA\tNA\t")
        end
        file.print("#{inter.type}\t#{inter.kegg}\t#{inter.url}\n")
    }
    STDERR.puts("wrote KEGGHPRD_result_reanalysis.txt")
}

# human_ligand
open(ARGV[0] + "/KEGGHPRD_result_cancer_ligand_reanalysis.txt", "w"){ |file|
    file.puts("Ligand\taverage count\treceptor ratio\tligand ratio\tactivate or inhibit\tpathway\tlink\treceptor(normalized count)")
    interactions.each{ |inter|
        file.print("#{inter.ligand_symbol}\t")

        sum_of_receptor_expression_human = 0.0
        sum_of_receptor_expression_mouse = 0.0
        interactions.each{ |inter2|
            if inter.ligand_symbol == inter2.ligand_symbol then
                sum_of_receptor_expression_human += inter2.receptor_human.expression if inter2.receptor_human != nil
                sum_of_receptor_expression_mouse += inter2.receptor_mouse.expression if inter2.receptor_mouse != nil
            end
        }
        if inter.ligand_human != nil then
            file.print("#{Math.sqrt((sum_of_receptor_expression_human + sum_of_receptor_expression_mouse) * inter.ligand_human.expression)}\t")
        else
            file.print("NA\t")
        end
        # receptor ratio
        if sum_of_receptor_expression_human + sum_of_receptor_expression_mouse > 0 then
            file.print("#{sum_of_receptor_expression_human / (sum_of_receptor_expression_human + sum_of_receptor_expression_mouse)}\t")
        else
            file.print("NA\t")
        end
        # ligand ratio
        if inter.ligand_human != nil && inter.ligand_mouse != nil && inter.ligand_human.expression + inter.ligand_mouse.expression > 0 then
            file.print("#{inter.ligand_human.expression / (inter.ligand_human.expression + inter.ligand_mouse.expression)}\t")
        else
            file.print("NA\t")
        end
        # type, kegg, url
        file.print("#{inter.type}\t#{inter.kegg}\t#{inter.url}\t")
        # receptor ranking
        receptors = []
        interactions.each{ |inter2|
            if inter.ligand_symbol == inter2.ligand_symbol then
                exp_sum = 0
                exp_sum += inter2.receptor_human.expression if inter2.receptor_human != nil
                exp_sum += inter2.receptor_mouse.expression if inter2.receptor_mouse != nil
                receptors << { receptor: inter2.receptor_symbol, exp: exp_sum }
            end
        }
        receptors.sort!{ |a, b|
            b[:exp] <=> a[:exp]
        }
        receptors.each{ |a|
            file.print("#{a[:receptor]}(#{a[:exp]})\t")
        }
        file.print("\n")
    }
    STDERR.puts("wrote KEGGHPRD_result_cancer_ligand_reanalysis.txt")
}

# mouse_ligand
open(ARGV[0] + "/KEGGHPRD_result_stroma_ligand_reanalysis.txt", "w"){ |file|
    file.puts("Ligand\taverage count\treceptor ratio\tligand ratio\tactivate or inhibit\tpathway\tlink\treceptor(normalized count)")
    interactions.each{ |inter|
        file.print("#{inter.ligand_symbol}\t")

        sum_of_receptor_expression_human = 0.0
        sum_of_receptor_expression_mouse = 0.0
        interactions.each{ |inter2|
            if inter.ligand_symbol == inter2.ligand_symbol then
                sum_of_receptor_expression_human += inter2.receptor_human.expression if inter2.receptor_human != nil
                sum_of_receptor_expression_mouse += inter2.receptor_mouse.expression if inter2.receptor_mouse != nil
            end
        }
        if inter.ligand_mouse != nil then
            file.print("#{Math.sqrt((sum_of_receptor_expression_human + sum_of_receptor_expression_mouse) * inter.ligand_mouse.expression)}\t")
        else
            file.print("NA\t")
        end
        # receptor ratio
        if sum_of_receptor_expression_human + sum_of_receptor_expression_mouse > 0 then
            file.print("#{sum_of_receptor_expression_mouse / (sum_of_receptor_expression_human + sum_of_receptor_expression_mouse)}\t")
        else
            file.print("NA\t")
        end
        # ligand ratio
        if inter.ligand_human != nil && inter.ligand_mouse != nil && inter.ligand_human.expression + inter.ligand_mouse.expression > 0 then
            file.print("#{inter.ligand_mouse.expression / (inter.ligand_human.expression + inter.ligand_mouse.expression)}\t")
        else
            file.print("NA\t")
        end
        # type, kegg, url
        file.print("#{inter.type}\t#{inter.kegg}\t#{inter.url}\t")
        # receptor ranking
        receptors = []
        interactions.each{ |inter2|
            if inter.ligand_symbol == inter2.ligand_symbol then
                exp_sum = 0
                exp_sum += inter2.receptor_human.expression if inter2.receptor_human != nil
                exp_sum += inter2.receptor_mouse.expression if inter2.receptor_mouse != nil
                receptors << { receptor: inter2.receptor_symbol, exp: exp_sum }
            end
        }
        receptors.sort!{ |a, b|
            b[:exp] <=> a[:exp]
        }
        receptors.each{ |a|
            file.print("#{a[:receptor]}(#{a[:exp]})\t")
        }
        file.print("\n")
    }
    STDERR.puts("wrote KEGGHPRD_result_stroma_ligand_reanalysis.txt")
}

# human_receptor
open(ARGV[0] + "/KEGGHPRD_result_cancer_receptor_reanalysis.txt", "w"){ |file|
    file.puts("Receptor\taverage count\tligand ratio\treceptor ratio\tactivate or inhibit\tpathway\tlink\tligand(normalized count)")
    interactions.each{ |inter|
        file.print("#{inter.receptor_symbol}\t")

        sum_of_ligand_expression_human = 0.0
        sum_of_ligand_expression_mouse = 0.0
        interactions.each{ |inter2|
            if inter.receptor_symbol == inter2.receptor_symbol then
                sum_of_ligand_expression_human += inter2.ligand_human.expression if inter2.ligand_human != nil
                sum_of_ligand_expression_mouse += inter2.ligand_mouse.expression if inter2.ligand_mouse != nil
            end
        }
        if inter.receptor_human != nil then
            file.print("#{Math.sqrt((sum_of_ligand_expression_human + sum_of_ligand_expression_mouse) * inter.receptor_human.expression)}\t")
        else
            file.print("NA\t")
        end
        # ligand ratio
        if sum_of_ligand_expression_human + sum_of_ligand_expression_mouse > 0 then
            file.print("#{sum_of_ligand_expression_human / (sum_of_ligand_expression_human + sum_of_ligand_expression_mouse)}\t")
        else
            file.print("NA\t")
        end
        # receptor ratio
        if inter.receptor_human != nil && inter.receptor_mouse != nil && inter.receptor_human.expression + inter.receptor_mouse.expression > 0 then
            file.print("#{inter.receptor_human.expression / (inter.receptor_human.expression + inter.receptor_mouse.expression)}\t")
        else
            file.print("NA\t")
        end
        # type, kegg, url
        file.print("#{inter.type}\t#{inter.kegg}\t#{inter.url}\t")
        # ligand ranking
        ligands = []
        interactions.each{ |inter2|
            if inter.receptor_symbol == inter2.receptor_symbol then
                exp_sum = 0
                exp_sum += inter2.ligand_human.expression if inter2.ligand_human != nil
                exp_sum += inter2.ligand_mouse.expression if inter2.ligand_mouse != nil
                ligands << { ligand: inter2.ligand_symbol, exp: exp_sum }
            end
        }
        ligands.sort!{ |a, b|
            b[:exp] <=> a[:exp]
        }
        ligands.each{ |a|
            file.print("#{a[:ligand]}(#{a[:exp]})\t")
        }
        file.print("\n")
    }
    STDERR.puts("wrote KEGGHPRD_result_cancer_receptor_reanalysis.txt")
}

# mouse_receptor
open(ARGV[0] + "/KEGGHPRD_result_stroma_receptor_reanalysis.txt", "w"){ |file|
    file.puts("Receptor\taverage count\tligand ratio\treceptor ratio\tactivate or inhibit\tpathway\tlink\tligand(normalized count)")
    interactions.each{ |inter|
        file.print("#{inter.receptor_symbol}\t")

        sum_of_ligand_expression_human = 0.0
        sum_of_ligand_expression_mouse = 0.0
        interactions.each{ |inter2|
            if inter.receptor_symbol == inter2.receptor_symbol then
                sum_of_ligand_expression_human += inter2.ligand_human.expression if inter2.ligand_human != nil
                sum_of_ligand_expression_mouse += inter2.ligand_mouse.expression if inter2.ligand_mouse != nil
            end
        }
        if inter.receptor_mouse != nil then
            file.print("#{Math.sqrt((sum_of_ligand_expression_human + sum_of_ligand_expression_mouse) * inter.receptor_mouse.expression)}\t")
        else
            file.print("NA\t")
        end
        # ligand ratio
        if sum_of_ligand_expression_human + sum_of_ligand_expression_mouse > 0 then
            file.print("#{sum_of_ligand_expression_mouse / (sum_of_ligand_expression_human + sum_of_ligand_expression_mouse)}\t")
        else
            file.print("NA\t")
        end
        # receptor ratio
        if inter.receptor_human != nil && inter.receptor_mouse != nil && inter.receptor_human.expression + inter.receptor_mouse.expression > 0 then
            file.print("#{inter.receptor_mouse.expression / (inter.receptor_human.expression + inter.receptor_mouse.expression)}\t")
        else
            file.print("NA\t")
        end
        # type, kegg, url
        file.print("#{inter.type}\t#{inter.kegg}\t#{inter.url}\t")
        # ligand ranking
        ligands = []
        interactions.each{ |inter2|
            if inter.receptor_symbol == inter2.receptor_symbol then
                exp_sum = 0
                exp_sum += inter2.ligand_human.expression if inter2.ligand_human != nil
                exp_sum += inter2.ligand_mouse.expression if inter2.ligand_mouse != nil
                ligands << { ligand: inter2.ligand_symbol, exp: exp_sum }
            end
        }
        ligands.sort!{ |a, b|
            b[:exp] <=> a[:exp]
        }
        ligands.each{ |a|
            file.print("#{a[:ligand]}(#{a[:exp]})\t")
        }
        file.print("\n")
    }
    STDERR.puts("wrote KEGGHPRD_result_stroma_receptor_reanalysis.txt")
}



include("../src/Bioinformatics.jl")
using Plots, Pkg, Test

@testset "Tests" begin

    @testset "alignments.jl" begin
        seq1 = Bioinformatics.Sequence("AND", Bioinformatics.AA)
        seq2 = Bioinformatics.Sequence("SEND", Bioinformatics.AA)
        sm = Bioinformatics.BLOSUM62
        mat = Bioinformatics.pairwise_global_alignment_linear_gap(
            seq1,
            seq2,
            sm,
            10
        )
        @test mat[length(seq1)+1, length(seq2)+1] == 3

        seq1 = Bioinformatics.Sequence("HEAGAWGHEE", Bioinformatics.AA)
        seq2 = Bioinformatics.Sequence("PAWHEAE", Bioinformatics.AA)
        sm = Bioinformatics.BLOSUM50
        mat = Bioinformatics.pairwise_global_alignment_linear_gap(
            seq1,
            seq2,
            sm,
            8
        )
        @test mat[length(seq1)+1, length(seq2)+1] == 1
    end

    @testset "distances.jl" begin
        s1 = "kitten"
        s2 = "sitting"
        @test Bioinformatics.edit_dist(s1, s2) == 3

        s3 = "karolin"
        s4 = "kathrin"
        @test Bioinformatics.hamming_dist(s3, s4) == 3
    end

    @testset "io.jl" begin
        seq = Bioinformatics.readFASTA(joinpath(
            Pkg.dir("Bioinformatics"),
            "example_data",
            "NC_000017.fasta"
        ))
        @test length(collect(seq)) == 1
    end

    @testset "plots.jl" begin
        ENV["PLOTS_TEST"] = "true"
        ENV["GKSwstype"] = "100"

        s1 = Bioinformatics.Sequence("CGATATAGATT", Bioinformatics.DNA)
        s2 = Bioinformatics.Sequence("TATATAGTAT", Bioinformatics.DNA)

        mat = Bioinformatics.dotmatrix(s1, s2)
        plot = Bioinformatics.plot_dotmatrix(mat)
        @test isa(plot, Plots.Plot) == true
        @test isa(display(plot), Nothing) == true

        seq = Bioinformatics.Sequence(
            "CATGGGCATCGGCCATACGCC",
            Bioinformatics.DNA
        )
        plot = Bioinformatics.skew_plot(seq)
        @test isa(plot, Plots.Plot) == true
        @test isa(display(plot), Nothing) == true

        plot = Bioinformatics.plot_gc_content(seq, 7)
        @test isa(plot, Plots.Plot) == true
        @test isa(display(plot), Nothing) == true
    end

    @testset "sequence.jl" begin
        seq = Bioinformatics.Sequence("ATGACAGAT", Bioinformatics.DNA)

        @test Bioinformatics.transcription(seq) == Bioinformatics.Sequence(
            "AUGACAGAU",
            Bioinformatics.RNA
        )
        @test Bioinformatics.reverse_complement(seq) == Bioinformatics.Sequence(
            "ATCTGTCAT",
            Bioinformatics.DNA
        )

        @test Bioinformatics.translation(seq) == Bioinformatics.Sequence(
            "MTD",
            Bioinformatics.AA
        )

        @test Bioinformatics.kmers(seq, 8) == ["ATGACAGA", "TGACAGAT"]

        @test collect(values(Bioinformatics.possible_proteins(seq)))[1] == Bioinformatics.Sequence(
            "MTD",
            Bioinformatics.AA
        )

        seq = Bioinformatics.Sequence("ACDEFGHIKLMNPQRSTVWY", Bioinformatics.AA)
        @test_throws ErrorException Bioinformatics.transcription(seq)
        @test_throws ErrorException Bioinformatics.reverse_complement(seq)
        @test_throws ErrorException Bioinformatics.translation(seq)
        @test_throws ErrorException Bioinformatics.reading_frames(seq)
        @test_throws ErrorException Bioinformatics.possible_proteins(seq)

    end

    @testset "stats.jl" begin
        seq = Bioinformatics.Sequence("atagataactcgcatag", Bioinformatics.DNA)

        freqs = Bioinformatics.frequency(seq)
        @test freqs['A'] == 7
        @test freqs['C'] == 3
        @test freqs['G'] == 3
        @test freqs['T'] == 4

        seq = Bioinformatics.Sequence(
            "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGG",
            Bioinformatics.DNA
        )
        gc_content = Bioinformatics.gc_content(seq)
        @test round(gc_content, digits = 2) == 0.61

        seq = Bioinformatics.Sequence(
            collect(values(Bioinformatics.readFASTA("../example_data/P35858.fasta")))[1],
            Bioinformatics.AA
        )
        seq_stats = Bioinformatics.protparam(seq)
        @test seq_stats["Number of amino acids"] == 605
        @test seq_stats["Molecular weight"] - 66035.0 < 0.01
        @test Bioinformatics.isoelectric_point(seq) - 6.35 < 0.01
        @test seq_stats["# of negatively charged residues"] == 59
        @test seq_stats["# of positively charged residues"] == 54
        @test seq_stats["Extinction coefficient"] == 61555
        @test seq_stats["Instability index"] - 46.09 < 0.01
        @test seq_stats["Aliphatic index"] - 111.652 < 0.01
        @test seq_stats["Grand average of hydropathicity (GRAVY)"] -
              0.018843 < 0.01
    end

end

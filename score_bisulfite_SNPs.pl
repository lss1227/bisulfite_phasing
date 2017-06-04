sub score_bisulfite_SNPs{
    my ($line,$chr,$genomic_pos,$read_base,$genome1,$genome2,$irrelevant) = @_;
    my $first_read_conversion;
    my $genome_conversion;
    # warn "$line\n"; sleep(1);

    my $strand;
    while ( $line =~ /(XR|XG):Z:([^\t]+)/g ) {
        my $tag = $1;
        my $value = $2;
        chomp $value;
        if ($tag eq "XR") {
            $first_read_conversion = $value;
            $first_read_conversion =~ s/\r//;
            }
        elsif ($tag eq "XG") {
            $genome_conversion = $value;
            $genome_conversion =~ s/\r//;
            }
            }
    

	if ($first_read_conversion eq 'CT' and $genome_conversion eq 'CT') {
        $strand = 'OT';		## this is OT
        }
    elsif ($first_read_conversion eq 'GA' and $genome_conversion eq 'CT') {
        $strand = 'CTOT';		## this is CTOT
        }
    elsif ($first_read_conversion eq 'GA' and $genome_conversion eq 'GA') {
        $strand = 'CTOB';		## this is CTOB
    }
    elsif ($first_read_conversion eq 'CT' and $genome_conversion eq 'GA') {
        $strand = 'OB';		## this is OB
    }
    else {
        die "Unexpected combination of read and genome conversion: $first_read_conversion / $genome_conversion\n";
    }

    ### if the SNP in question involves a genomic C position we potentially need to allow both C or T as C (G/A for the reverse strand)

    if ( $snps{$chr}->{$genomic_pos}->{ref} eq 'C'){    # the position is a C in the reference genome
        ### C>T SNP. Can't use the position for top strand alignments
        if ($snps{$chr}->{$genomic_pos}->{snp} eq 'T'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($verbose){ warn "This is a C>T SNP and can't be considered for OT or CTOT alignments!\n";}
	            ++$ct_snp;
                return ($genome1,$genome2,$irrelevant);
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'C'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'T'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
                }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'G'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'C' or $read_base eq 'T'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'G'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'C'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'G' or $read_base eq 'A'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'A'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'C' or $read_base eq 'T'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'A'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'C'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'A'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'C') {
            die "The ref is the same as snp, which should not happen in snpfile.\n"
        }
    }

    if ( $snps{$chr}->{$genomic_pos}->{ref} eq 'A'){ 
        if ($snps{$chr}->{$genomic_pos}->{snp} eq 'G'){
            if ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($verbose){ warn "This is a A>G SNP and can't be considered for OB or CTOB alignments!\n";}
	            ++$ct_snp;
                return ($genome1,$genome2,$irrelevant);
            }
            elsif ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'A'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'G'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'C'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'A'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'C' or $read_base eq 'T'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'A'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'C'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'T'){
            if ($read_base eq 'A'){
                if ($verbose){ warn "genome 1 specific!\n";}
                ++$genome1;
            }
            elsif ($read_base eq 'T'){
                if ($verbose){ warn "genome 2 specific!\n";}
                ++$genome2;
            }
            else{
                if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
                ++$irrelevant;
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'A') {
            die "The ref is the same as snp, which should not happen in snpfile.\n"
        }
    }

    if ( $snps{$chr}->{$genomic_pos}->{ref} eq 'T'){ 
        if ($snps{$chr}->{$genomic_pos}->{snp} eq 'C'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($verbose){ warn "This is a T>C SNP and can't be considered for OT or CTOT alignments!\n";}
	            ++$ct_snp;
                return ($genome1,$genome2,$irrelevant);
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'T'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'C'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'G'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'T'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'G'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'T'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'G' or $read_base eq 'A'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'A'){
            if ($read_base eq 'T'){
                if ($verbose){ warn "genome 1 specific!\n";}
                ++$genome1;
            }
            elsif ($read_base eq 'A'){
                if ($verbose){ warn "genome 2 specific!\n";}
                ++$genome2;
            }
            else{
                if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
                ++$irrelevant;
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'T') {
            die "The ref is the same as snp, which should not happen in snpfile.\n"
        }
    }
        
    if ( $snps{$chr}->{$genomic_pos}->{ref} eq 'G'){
        if ($snps{$chr}->{$genomic_pos}->{snp} eq 'A'){
            if ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($verbose){ warn "This is a G>A SNP and can't be considered for OB or CTOB alignments!\n";}
	            ++$ct_snp;
                return ($genome1,$genome2,$irrelevant);
            }
            elsif ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'G'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'T'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'T'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'G'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'T'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'G' or $read_base eq 'A'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'T'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'C'){
            if ($strand eq 'OT' or $strand eq 'CTOT'){
                if ($read_base eq 'C' or $read_base eq 'T'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                elsif ($read_base eq 'G'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            elsif ($strand eq 'OB' or $strand eq 'CTOB'){
                if ($read_base eq 'G' or $read_base eq 'A'){
                    if ($verbose){ warn "genome 1 specific!\n";}
                    ++$genome1;
                }
                elsif ($read_base eq 'C'){
                    if ($verbose){ warn "genome 2 specific!\n";}
                    ++$genome2;
                }
                else{
                    if ($verbose){ warn "The base was different than both genome 1 or genome 2: $read_base\n";}
	                ++$irrelevant;
                }
            }
            else{
                die "Unable to identify read alignment strand\n";
            }
        }
        elsif ($snps{$chr}->{$genomic_pos}->{snp} eq 'G') {
            die "The ref is the same as snp, which should not happen in snpfile.\n"
        }
    }

    return ($genome1,$genome2,$irrelevant);
    
}
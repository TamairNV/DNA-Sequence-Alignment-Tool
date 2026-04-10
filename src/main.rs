use rand::prelude::*;
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::time::Instant;
use needletail::parse_fastx_file;
use std::fs::File;
use std::io::{Write, BufWriter};

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
struct Sequence{
    data : Vec<u64>,
    pub len : usize
}

impl Sequence {

    #[inline]
    pub fn get(&self, index: usize) -> u8 {
        let block_idx = index / 32;
        let bit_offset = (index % 32) * 2;
        // Shift the specific 2 bits to the far right, then mask with 0b11 (3)
        ((self.data[block_idx] >> bit_offset) & 0b11) as u8
    }

    pub fn create_random(&mut self) {
        let mut rng = rand::thread_rng();

        for block in self.data.iter_mut() {

            *block = rng.r#gen();
        }
        println!("Created random DNA");
    }


    pub fn from_fasta(file_path: &str) -> Self {
        let mut reader = parse_fastx_file(file_path).expect("Failed to find or open FASTA file");

        let mut packed_data: Vec<u64> = Vec::new();
        let mut total_len: usize = 0;

        let mut current_block: u64 = 0;
        let mut bases_in_block = 0;

        // Loop through the sequences in the file
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid FASTA record");

            for &base in seqrec.seq().iter() {


                let bit_val: u64 = match base {
                    b'A' | b'a' => 0b00,
                    b'C' | b'c' => 0b01,
                    b'G' | b'g' => 0b10,
                    b'T' | b't' => 0b11,
                    _ => 0b00, // Catch unknown bases and just default them to A
                };


                current_block |= bit_val << (bases_in_block * 2);
                bases_in_block += 1;
                total_len += 1;


                if bases_in_block == 32 {
                    packed_data.push(current_block);
                    current_block = 0;
                    bases_in_block = 0;
                }
            }
        }

        if bases_in_block > 0 {
            packed_data.push(current_block);
        }

        println!("Successfully loaded {} bases from FASTA!", total_len);

        Self {

            data: packed_data,
            len:total_len,
        }



    }



    pub fn create_table(&self, kmer: usize,w:usize) -> Vec<(u64, usize,bool)> {
        let mut kmer_locations: Vec<(u64, usize,bool)> = Vec::new();

        let mut window_buffer: VecDeque<(u64, usize)> = VecDeque::with_capacity(4);

        let mut last_saved_pos = usize::MAX;

        let bits_per_kmer = kmer * 2;

        //  overflow check
        let mask: u64 = if bits_per_kmer >= 64 {
            !0 // This is all 1s (0xFFFFFFFFFFFFFFFF)
        } else {
            (1 << bits_per_kmer) - 1
        };

        let mut window: u64 = 0;

        for i in 0..self.len {
            let new_piece = self.get(i) as u64;

            //sliding window
            window = (window << 2) | new_piece;
            window &= mask;

            if i >= kmer - 1 {
                let start_pos = i - (kmer - 1);
                window_buffer.push_back((window, start_pos));

                if window_buffer.len() > w {
                    window_buffer.pop_front();
                }

                if window_buffer.len() == w {
                    let min_item = window_buffer.iter().min_by_key(|&(k, _)| k).unwrap();
                    let (min_kmer, min_pos) = *min_item;

                    if min_pos != last_saved_pos {

                        let rev = reverse_complement(min_kmer, kmer);

                        if rev < min_kmer {
                            kmer_locations.push((rev, min_pos,true));
                        }else{
                            kmer_locations.push((min_kmer, min_pos,false));
                        }



                        last_saved_pos = min_pos;
                    }
                }
            }
        }
        kmer_locations
    }



    pub fn to_string(&self)->String{
        let mut str = String::from("");
        for &block in self.data.iter() {
            for i in 0..32{
                let shift = i * 2;
                let item = (block >> shift) & 0b11;

                if item == 0{
                    str.push('A');
                }
                else if item ==1{
                    str.push('T');
                }
                else if item ==2{
                    str.push('C');
                }
                else if item ==3{
                    str.push('G');
                }

            }
        }
        str
    }


}
pub fn reverse_complement(kmer: u64, k: usize) -> u64 {
    let mut rc = 0;
    for i in 0..k {

        let base = (kmer >> (i * 2)) & 0b11;

        let complement = base ^ 0b11;

        rc = (rc << 2) | complement;
    }
    rc
}

pub fn two_pointer_matcher(dna_array_a: &[(u64, usize,bool)], dna_array_b: &[(u64, usize,bool)]) -> Vec<(usize, usize, bool)> {
    let mut pointer_a = 0;
    let mut pointer_b = 0;
    let mut final_seeds = Vec::new();

    //. If a k-mer appears more than this, it's skipped.
    let max_frequency = 5;

    while pointer_a < dna_array_a.len() && pointer_b < dna_array_b.len() {
        let kmer_a = dna_array_a[pointer_a].0;
        let kmer_b = dna_array_b[pointer_b].0;

        if kmer_a < kmer_b {
            pointer_a += 1;
        } else if kmer_b < kmer_a {
            pointer_b += 1;
        } else {
            let matched_kmer = kmer_a;

            let mut end_a = pointer_a;
            while end_a < dna_array_a.len() && dna_array_a[end_a].0 == matched_kmer {
                end_a += 1;
            }

            let mut end_b = pointer_b;
            while end_b < dna_array_b.len() && dna_array_b[end_b].0 == matched_kmer {
                end_b += 1;
            }

            let count_a = end_a - pointer_a;
            let count_b = end_b - pointer_b;

            if count_a <= max_frequency && count_b <= max_frequency {
                for i in pointer_a..end_a {
                    for j in pointer_b..end_b {
                        let is_reversed_match = dna_array_a[i].2 ^ dna_array_b[j].2; // XOR
                        final_seeds.push((dna_array_a[i].1, dna_array_b[j].1, is_reversed_match));
                    }
                }
            }

            pointer_a = end_a;
            pointer_b = end_b;
        }
    }

    final_seeds
}
#[derive(Debug,Clone)]
pub enum CigarOp {
    Match(u16),  // 'M'
    Insert(u16), // 'I'
    Delete(u16), // 'D'
}

#[derive(Debug,Clone)]
pub struct AlignmentResult {
    pub score: i16,
    pub a_start: usize,
    pub a_end: usize,
    pub b_start: usize,
    pub b_end: usize,
    pub cigar: Vec<CigarOp>,
}

fn smith_waterman(seq_a : &Vec<u8>, seq_b : &Vec<u8>, max_score : i16, match_score:i16, mismatch_score:i16, gap_score:i16, start_a : usize, start_b : usize) -> (usize, usize, Vec<CigarOp>) {


    let mut dna_start_a = start_a;
    let mut dna_start_b = start_b;
    let mut dna_end_a = 0;
    let mut dna_end_b = 0;

    let mut r_index = 1;
    let mut c_index = 1;

    let row_length = seq_a.len();
    let col_length = seq_b.len() ;
    let table_size = (row_length + 1) * (col_length + 1);
    let iter_amount = (row_length) * (col_length);
    let mut table: Vec<i16> = vec![0; table_size];
    let mut current_position = 0;
    for i in 0..iter_amount {

        if r_index >= row_length+1 {
            r_index = 1;
            c_index +=1;
        }

        let r_index_2d = r_index + c_index * (row_length + 1);

        let left = &table[r_index_2d -1]  + gap_score ;

        let top = &table[r_index_2d - (row_length + 1)] + gap_score;

        let mut current_match_score = mismatch_score;

        if seq_a[r_index-1] == seq_b[c_index-1]{
            current_match_score = match_score;
        }
        let diag = &table[r_index_2d - (row_length + 1) - 1] + current_match_score as i16;
        let score = 0.max(diag).max(top.max(left));
        table[r_index_2d] = score;
        if score == max_score {
            current_position = r_index_2d;
            break;
        }

        r_index += 1;
    }

    let mut cigar =  Vec::<CigarOp>::new();
    let mut current_type : u8 = 0;
    let mut previous_type:u8  = 0;
    let mut count = 1;
    let width = row_length + 1;

    let mut r_index = current_position % width;
    let mut c_index = current_position / width;

    dna_end_a = start_a + r_index;
    dna_end_b = start_b + c_index;

    while table[current_position] > 0{
        if current_position == 0{
            break;
        }
        if r_index == 0 || c_index == 0 {
            break;
        }
        r_index = current_position % width;
        c_index = current_position / width;

        let is_match = seq_a[r_index - 1] == seq_b[c_index - 1];

        let mut diag_points : i16 = mismatch_score;
        if is_match {
            diag_points = match_score;
        }

        let diag_index = current_position - width - 1;
        let diag = &table[diag_index] + diag_points;
        if diag == table[current_position] {
            current_position = diag_index;
            current_type = 1;

        }else{
            let left_index = current_position - 1;
            let left = &table[left_index] + gap_score;
            if left == table[current_position] {
                current_position = left_index;
                current_type = 2;

            }else {
                let top_index = current_position - width;
                let top = &table[top_index] + gap_score;
                if top == table[current_position] {
                    current_position = top_index;
                    current_type = 3;

                }
            }
        }



        if previous_type == 0 {
            previous_type = current_type;
            count = 1;
        } else if current_type == previous_type {

            count += 1;
        } else {

            match previous_type {
                1 => cigar.push(CigarOp::Match(count)),
                2 => cigar.push(CigarOp::Insert(count)),
                3 => cigar.push(CigarOp::Delete(count)),
                _ => unreachable!(), // 0 should never reach here
            }


            previous_type = current_type;
            count = 1;
        }





    }
    if count > 0 {
        match previous_type {
            1 => cigar.push(CigarOp::Match(count)),
            2 => cigar.push(CigarOp::Insert(count)),
            3 => cigar.push(CigarOp::Delete(count)),
            _ => {}
        }
    }

    dna_start_a = start_a + r_index;
    dna_start_b = start_b + c_index;

    cigar.reverse();

    let result = (dna_end_a,dna_end_b,cigar);


    result

}

pub fn two_row_smith_waterman_chunked<'a>(
    seq_a: &[u8],
    seq_b: &[u8],
    mut top_row: &'a mut [i16],
    mut bottom_row: &'a mut [i16],
    match_score: i16,
    mismatch_score: i16,
    gap_score: i16,
    mut current_max: i16
) -> (i16, bool){

    for i in 1..=seq_b.len() {

        bottom_row[0] = 0;

        for j in 1..=seq_a.len() {

            let left = bottom_row[j - 1] + gap_score;
            let top = top_row[j] + gap_score;

            let current_match_score = if seq_a[j - 1] == seq_b[i - 1] {
                match_score
            } else {
                mismatch_score
            };

            let diag = top_row[j - 1] + current_match_score;

            let score = 0.max(diag).max(top.max(left));

            bottom_row[j] = score;

            if score > current_max {
                current_max = score;
            }
        }

        std::mem::swap(&mut top_row, &mut bottom_row);
    }
    let hit_edge = bottom_row[seq_a.len()] > 0;

    (current_max, hit_edge)

}
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
pub fn expand_seeds(seeds : &Vec<(usize, usize, bool)>, dna_a : &Sequence, dna_b : &Sequence, chunk_size  : usize, kmer: usize, match_score: i16,
                    mismatch_score: i16,
                    gap_score: i16,) -> Vec<AlignmentResult> {

    let pb = ProgressBar::new(seeds.len() as u64);

    // 2. Set a custom style
    pb.set_style(
        ProgressStyle::with_template(
            // {wide_bar} automatically adjusts to your terminal width
            // {eta} gives you the estimated time remaining
            "Expanding seeds: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} (ETA: {eta})"
        )
            .unwrap()
            // These characters dictate how the bar looks: "filled > empty"
            .progress_chars("=>-")
    );

    let results: Vec<AlignmentResult> = seeds.par_iter()
        .progress_with(pb) // <-- Changed this!
        .map(|seed| {
            expand_seed(seed, dna_a, dna_b, chunk_size, match_score, mismatch_score, gap_score, kmer)
        })
        .collect();
    results

}

fn expand_seed(seed : &(usize, usize, bool), dna_a : &Sequence, dna_b : &Sequence, chunk_size  : usize, match_score: i16,
               mismatch_score: i16,
               gap_score: i16,
               kmer : usize) -> AlignmentResult{
    let mut top_row = vec![0; chunk_size + 1];
    let mut bottom_row = vec![0; chunk_size + 1];

    let max_extension = 5000; // Cap at 3000 bases
    let mut right_a_pos = seed.0;
    let mut right_b_pos = seed.1;
    let mut right_score = kmer as i16 * match_score;

    loop {

        let end_a = (right_a_pos + chunk_size).min(dna_a.len);
        let end_a = (right_a_pos + chunk_size).min(dna_a.len);

        //if reversed
        let (start_b, end_b, chunk_b) = if seed.2 {
            let start = right_b_pos.saturating_sub(chunk_size);
            let end = right_b_pos;
            let mut b_vec = unpack_region(dna_b, start, end);

            b_vec = manually_reverse_complement_vec(&b_vec);
            (start, end, b_vec)
        } else {

            let start = right_b_pos;
            let end = (right_b_pos + chunk_size).min(dna_b.len);
            let b_vec = unpack_region(dna_b, start, end);
            (start, end, b_vec)
        };

        let chunk_a = unpack_region(dna_a, right_a_pos, end_a);

        let (new_max, hit_edge) = two_row_smith_waterman_chunked(
            &chunk_a, &chunk_b, &mut top_row, &mut bottom_row,
            match_score, mismatch_score, gap_score, right_score,
        );
        right_score = new_max;
        let current_distance = right_a_pos.saturating_sub(seed.0);
        if hit_edge && current_distance < max_extension {
            right_a_pos += chunk_size;

            if seed.2 {
                right_b_pos = start_b; // move left in B
                if start_b == 0 { break; } // Hit the edge
            } else {
                right_b_pos = end_b; // move right in B
                if end_b >= dna_b.len { break; }
            }

            if end_a >= dna_a.len { break; }
        } else {
            break;
        }
    }


    let final_right_end_a = (right_a_pos + chunk_size).min(dna_a.len);

    let right_seq_a = unpack_region(dna_a, seed.0 + kmer, final_right_end_a);
    let (final_right_b_start, final_right_b_end) = if seed.2 {
        (right_b_pos, seed.1)
    } else {
        (seed.1 + kmer, (right_b_pos + chunk_size).min(dna_b.len))
    };

    let right_seq_b = unpack_region(dna_b, final_right_b_start, final_right_b_end);


    let right_result = smith_waterman(&right_seq_a, &right_seq_b, right_score, match_score, mismatch_score, gap_score, seed.0, seed.1);

    top_row.fill(0);
    bottom_row.fill(0);

    let mut left_a_pos = seed.0;
    let mut left_b_pos = seed.1;
    let mut left_score = kmer as i16 * match_score;
    loop {
        let start_a = left_a_pos.saturating_sub(chunk_size);
        let mut chunk_a = unpack_region(dna_a, start_a, left_a_pos);
        chunk_a.reverse(); // Standard "walking backward" trick for A

        //if reversed inverse B
        let (start_b, end_b, chunk_b) = if seed.2 {

            let start = left_b_pos;
            let end = (left_b_pos + chunk_size).min(dna_b.len);
            let b_vec = unpack_region(dna_b, start, end);


            let chunk_b_rc = manually_reverse_complement_vec(&b_vec);
            (start, end, chunk_b_rc)
        } else {

            let start = left_b_pos.saturating_sub(chunk_size);
            let end = left_b_pos;
            let mut b_vec = unpack_region(dna_b, start, end);
            b_vec.reverse();
            (start, end, b_vec)
        };

        let (new_max, hit_edge) = two_row_smith_waterman_chunked(
            &chunk_a, &chunk_b, &mut top_row, &mut bottom_row,
            match_score, mismatch_score, gap_score, left_score,
        );

        left_score = new_max;
        let current_left_distance = seed.0.saturating_sub(left_a_pos);

        if hit_edge && current_left_distance < max_extension {

            left_a_pos = start_a;
            if start_a == 0 { break; }


            if seed.2 {
                left_b_pos = end_b;
                if end_b >= dna_b.len { break; }
            } else {
                left_b_pos = start_b;
                if start_b == 0 { break; }
            }
        } else {
            break;
        }
    }

    let left_start_a = left_a_pos.saturating_sub(chunk_size);

    let mut left_seq_a = unpack_region(dna_a, left_start_a, seed.0);
    let (final_left_b_start, final_left_b_end) = if seed.2 {
        (seed.1 + kmer, left_b_pos)
    } else {
        (left_b_pos.saturating_sub(chunk_size), seed.1)
    };

    left_seq_a.reverse();
    left_seq_a.reverse();
    // Pass 0, 0 so the result is just the raw length of the backward match
    let mut left_result = smith_waterman(&left_seq_a, &left_seq_a, left_score, match_score, mismatch_score, gap_score, 0, 0);
    left_result.2.reverse();

    let mut final_cigar = Vec::new();
    final_cigar.extend(left_result.2);
    final_cigar.push(CigarOp::Match(kmer as u16));
    final_cigar.extend(right_result.2);
    let (final_b_start, final_b_end) = if seed.2 {
        let b_start = right_result.1;
        let b_end = left_result.1;

        (b_start, b_end)
    } else {

        let b_start = seed.1.saturating_sub(left_result.1);
        let b_end = right_result.1;
        (b_start, b_end)
    };

    let result = AlignmentResult {
        score: right_score + left_score - (kmer as i16 * match_score),
        // A is always the forward reference
        a_start: seed.0.saturating_sub(left_result.0),
        a_end: right_result.0,

        b_start: final_b_start,
        b_end: final_b_end,

        cigar: final_cigar,
    };

    result

}

pub fn manually_reverse_complement_vec(chunk: &[u8]) -> Vec<u8> {
    chunk.iter()
        .rev()
        .map(|&ch| ch ^ 0b11)
        .collect()
}

fn unpack_region(seq: &Sequence, start: usize, end: usize) -> Vec<u8> {
    if start >= end {
        return Vec::new();
    }

    if start >= seq.len {
        return Vec::new();
    }

    let mut unpacked = Vec::with_capacity(end - start);
    for i in start..end {
        unpacked.push(seq.get(i));
    }
    unpacked
}



use rustc_hash::FxHashMap;

pub fn remove_duplicates(mut alignments: Vec<AlignmentResult>) -> Vec<AlignmentResult> {
    alignments.sort_unstable_by(|x, y| y.score.cmp(&x.score));

    let mut final_results: Vec<AlignmentResult> = Vec::new();

    let bucket_size = 50_000;

    let mut grid: FxHashMap<usize, Vec<usize>> = FxHashMap::default();

    for current in alignments {
        let mut is_duplicate = false;

        let start_bucket = current.a_start / bucket_size;
        let end_bucket = current.a_end / bucket_size;

        'check: for bucket in start_bucket..=end_bucket {
            if let Some(kept_indices) = grid.get(&bucket) {
                for &idx in kept_indices {
                    let kept = &final_results[idx];
                    let overlap_a = current.a_start < kept.a_end && current.a_end > kept.a_start;
                    let overlap_b = current.b_start < kept.b_end && current.b_end > kept.b_start;

                    if overlap_a && overlap_b {
                        is_duplicate = true;
                        break 'check;
                    }
                }
            }
        }

        if !is_duplicate {
            let new_idx = final_results.len();
            for bucket in start_bucket..=end_bucket {
                grid.entry(bucket).or_default().push(new_idx);
            }
            final_results.push(current);
        }
    }

    final_results
}

pub fn alignments_to_csv(alignments: &[AlignmentResult], file_path: &str) -> std::io::Result<()> {
    let file = File::create(file_path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "a_start,a_end,b_start,b_end,score,length,cigar")?;

    for reg in alignments {
        // Only export significant matches to keep the file size down
        if reg.score < 100 { continue; }

        let length = reg.a_end - reg.a_start;
        let cigar_str = format_cigar(&reg.cigar);

        writeln!(
            writer,
            "{},{},{},{},{},{},{}",
            reg.a_start,
            reg.a_end,
            reg.b_start,
            reg.b_end,
            reg.score,
            length,
            cigar_str
        )?;
    }

    writer.flush()?;
    println!("CSV exported to: {}", file_path);
    Ok(())
}

fn format_cigar(ops: &[CigarOp]) -> String {
    let mut s = String::new();
    for op in ops {
        match op {
            CigarOp::Match(len) => s.push_str(&format!("{}M", len)),
            CigarOp::Insert(len) => s.push_str(&format!("{}I", len)),
            CigarOp::Delete(len) => s.push_str(&format!("{}D", len)),
        }
    }
    s
}
pub fn get_alignments(k_mer: usize, w:usize, match_score: i16,
                      mismatch_score: i16,
                      gap_score: i16, genome1: Box<str>,genome2: Box<str>) -> Vec<AlignmentResult> {

    let mut dna_a = Sequence::from_fasta(&genome1);
    let mut dna_b = Sequence::from_fasta(&*genome2);

    let start = Instant::now();

    let (mut dna_a_table, mut dna_b_table) = rayon::join(
        || dna_a.create_table(k_mer,w),
        || dna_b.create_table(k_mer,w)
    );
    println!("Created seed arrays");

    rayon::join(
        || dna_a_table.par_sort_unstable(),
        || dna_b_table.par_sort_unstable()
    );
    println!("sorted kmer arrays" );
    let comparison_locations = two_pointer_matcher(&dna_a_table,&dna_b_table);
    println!("created comparison array");


    println!("expanding seeds");

    let mut alignments = expand_seeds(&comparison_locations, &dna_a, &dna_b, 100, k_mer, match_score, mismatch_score, gap_score);
    println!("created full alignments");

    alignments.retain(|x| x.a_end.saturating_sub(x.a_start) > 800);

    let result = remove_duplicates(alignments);
    println!("removed dupes");
    println!("Done. Time elapsed: {:?}", start.elapsed());
    result


}
fn main() {
    let table = get_alignments(21, 50, 2, -1, -1,
                               Box::from("turnip/ncbi_dataset/data/GCF_000695525.1/GCF_000695525.1_BOL_genomic.fna"),
                               Box::from("thaliana/ncbi_dataset/data/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna"));
    alignments_to_csv(&*table, "data.csv");
    //println!("Comparison table score: {:?}", table);
}
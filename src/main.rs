use rand::prelude::*;
use rayon::prelude::*;
use rustc_hash::FxHashMap;


use std::collections::{HashMap, VecDeque};
struct Sequence{
    data : Vec<u64>,
    pub len : usize
}


impl Sequence {
    pub fn new(len: usize) -> Self {
        let num_blocks = (len + 31) / 32;
        Self {

            data: vec![0u64; num_blocks],
            len,
        }
    }
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



    pub fn create_table(&self, kmer: usize) -> Vec<(u64, usize)> {
        let mut kmer_locations: Vec<(u64, usize)> = Vec::new();

        // VecDeque is faster for "sliding" than a normal Vec
        let mut window_buffer: VecDeque<(u64, usize)> = VecDeque::with_capacity(4);
        let w = 10;
        let mut last_saved_pos = usize::MAX;

        let bits_per_kmer = kmer * 2;

        // Safety check for the 64-bit overflow
        let mask: u64 = if bits_per_kmer >= 64 {
            !0 // This is shorthand for "all 1s" (0xFFFFFFFFFFFFFFFF)
        } else {
            (1 << bits_per_kmer) - 1
        };

        let mut window: u64 = 0;

        for i in 0..self.len {
            let new_piece = self.get(i) as u64;

            // Standard sliding logic
            window = (window << 2) | new_piece;
            window &= mask;

            if i >= kmer - 1 {
                let start_pos = i - (kmer - 1);
                window_buffer.push_back((window, start_pos));

                if window_buffer.len() > w {
                    window_buffer.pop_front(); // Instant O(1) removal
                }

                if window_buffer.len() == w {
                    let min_item = window_buffer.iter().min_by_key(|&(k, _)| k).unwrap();
                    let (min_kmer, min_pos) = *min_item;

                    if min_pos != last_saved_pos {
                        kmer_locations.push((min_kmer, min_pos));

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

pub fn two_pointer_matcher(dna_array_a: &[(u64, usize)], dna_array_b: &[(u64, usize)]) -> Vec<(usize, usize)> {
    let mut pointer_a = 0;
    let mut pointer_b = 0;
    let mut final_seeds = Vec::new();

    // Set a limit. If a k-mer appears more than this, it's junk. Skip it.
    let max_frequency = 50;

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

            // THE FIX: Check how big the "clump" is before running the nested loop
            let count_a = end_a - pointer_a;
            let count_b = end_b - pointer_b;

            if count_a <= max_frequency && count_b <= max_frequency {
                for i in pointer_a..end_a {
                    for j in pointer_b..end_b {
                        final_seeds.push((dna_array_a[i].1, dna_array_b[j].1));
                    }
                }
            }

            pointer_a = end_a;
            pointer_b = end_b;
        }
    }
    final_seeds
}


fn smith_waterman(seq_a : &Vec<u8>, seq_b : &Vec<u8>,match_score:i16,mismatch_score:i16,gap_score:i16) -> (Vec<u8>,Vec<u8>) {

    let max_score = two_row_smith_waterman(&seq_a, &seq_b, match_score, mismatch_score, gap_score);
    let mut new_a_seq = vec![];
    let mut new_b_seq = vec![];
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

        let left = &table[r_index_2d -1]  + (gap_score * (c_index) as i16);

        let top = &table[r_index_2d - (row_length + 1)] + (gap_score * (r_index) as i16);

        let mut current_match_score = mismatch_score;

        if seq_a[r_index-1] == seq_b[c_index-1]{
            current_match_score = match_score;
        }
        let diag = &table[r_index_2d - (row_length + 1) - 1] + current_match_score as i16;
        let score = 0.max(diag).max(top.max(left));
        table[r_index_2d] = score;

        //println!("left {} + {} = {}",&table[r_index_2d -1],(gap_score * (c_index) as i16),left);
        //println!("top {} + {} = {}",&table[r_index_2d - (row_length + 1)],(gap_score * (r_index) as i16),top );
        //println!("diag {} + {} = {}",&table[r_index_2d - (row_length + 1) - 1],current_match_score,diag);
        //println!("-------");
        if score == max_score {
            current_position = r_index_2d;
            break;
        }

        r_index += 1;
    }


    while current_position != 0{
        if current_position == 0{
            break;
        }


        let width = row_length + 1;

        let r_index = current_position % width;
        let c_index = current_position / width;

        let is_match = seq_a[r_index - 1] == seq_b[c_index - 1];

        let mut diag_points : i16 = mismatch_score;
        if is_match {
            diag_points = match_score;
        }

        let diag_index = current_position - width - 1;
        let diag = &table[diag_index] + diag_points;
        if diag == table[current_position] {
            current_position = diag_index;
            new_a_seq.push(seq_a[r_index - 1]);
            new_b_seq.push(seq_b[c_index - 1]);

            continue;
        }

        let left_index = current_position - 1;
        let left = &table[left_index] + (gap_score * c_index as i16);
        if left == table[current_position] {
            current_position = left_index;
            //addition on seq b so dash on a
            new_b_seq.push(5);
            new_a_seq.push(seq_a[r_index - 1]);
            continue;
        }

        let top_index = current_position - width;
        let top = &table[top_index] + (gap_score * r_index as i16);
        if top == table[current_position] {
            current_position = top_index;
            //addition on seq a so dash on b
            new_a_seq.push(5);
            new_b_seq.push(seq_b[c_index - 1]);
            continue;
        }



    }

    new_a_seq.reverse();
    new_b_seq.reverse();
    (new_a_seq,new_b_seq)

}

fn two_row_smith_waterman(seq_a : &Vec<u8>, seq_b : &Vec<u8>,match_score:i16,mismatch_score:i16,gap_score:i16) -> i16 {

    let mut r_index = 1;
    let mut c_index = 1;

    let mut top_row = vec![0; seq_a.len() + 1];

    let mut bottom_row = vec![0; seq_a.len() + 1];
    let mut max_score = 0;
    let iter_amount = (seq_a.len()) * (seq_b.len());

    for _ in 0..iter_amount {

        if r_index >= seq_a.len()+1 {
            r_index = 1;
            c_index +=1;
            std::mem::swap(&mut top_row, &mut bottom_row);
            bottom_row.fill(0);
        }

        let left = &bottom_row[r_index-1]  + (gap_score * (c_index) as i16);

        let top = &top_row[r_index]+  (gap_score * (r_index) as i16);
        let mut current_match_score = mismatch_score;

        if seq_a[r_index-1] == seq_b[c_index-1]{
            current_match_score = match_score;
        }
        let diag = &top_row[r_index-1] + current_match_score;
        let score = 0.max(diag).max(top.max(left));

        bottom_row[r_index] = score;
        if score > max_score{
            max_score = score;
        }
        r_index += 1;

    }
    println!("{}",max_score);
    max_score
}
pub fn create_tables() -> Vec<(usize, usize)> {

    let mut dna_a = Sequence::new(1000000000);
    let mut dna_b = Sequence::new(1000000000);
    dna_a.create_random();
    dna_b.create_random();
    let start = Instant::now();
    let k_mer : usize = 21;
    let (mut dna_a_table, mut dna_b_table) = rayon::join(
        || dna_a.create_table(k_mer),
        || dna_b.create_table(k_mer)
    );
    println!("Created kmer arrays");

    rayon::join(
        || dna_a_table.par_sort_unstable(),
        || dna_b_table.par_sort_unstable()
    );
    println!("sorted kmer arrays");
    let result = two_pointer_matcher(&dna_a_table,&dna_b_table);
    println!("created comparison table");
    println!("Done. Time elapsed: {:?}", start.elapsed());
    result


}
use std::time::Instant;
use rayon::scope;

fn main() {

    //create_tables();

    let seq_a = vec![2,3,1,4,2,3,3,3,1,3,4,2,3,3,1,3,2];
    let seq_b = vec![2,3,1,4,2,3,3,1,4,2,3,1,3,2];


    let seq = smith_waterman(&seq_a,&seq_b,2,-1,-1);
    println!("Comparison table score: {:?} {:?}", seq.0, seq.1);


}

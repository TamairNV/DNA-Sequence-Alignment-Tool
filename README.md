# DNA Aligner from Scratch

A high-performance DNA sequence aligner built in Rust to explore the intersection of biology and big data. This project focuses on efficient memory management and algorithmic optimization for processing genomic data.

##  The Mission
Real value lies in the niche. This project was built to bridge the gap between computer science and bioinformatics, moving from a Python prototype to a high-speed Rust implementation in just two days.

##  Technical Deep Dive

### 1. Seed and Extend Strategy
Aligning massive genomes is a search problem. Instead of a brute-force comparison, this tool uses:
* **Seeding (Minimizers):** DNA is broken into short chunks. Minimizers allow for compressed indexing without losing accuracy.
* **Bit-Packing:** DNA bases are represented in just **2 bits**, packing 32-base chunks into a single 64-bit integer. This enables binary operator processing in just a few clock cycles.

### 2. Optimized Smith-Waterman Alignment
Once a matching seed is found, it expands using the Smith-Waterman algorithm with a custom scoring matrix (matching rewards vs. gap penalties).

### 3. The Two-Pass Memory Hack
Standard alignment requires an $O(N \times M)$ grid, which is incredibly memory-heavy. This implementation uses a custom two-pass system to keep the footprint tiny:
* **Score Pass:** Uses a sliding two-row memory technique. By only keeping the current and previous rows in RAM, the algorithm identifies the "peak" score coordinate without storing the entire table.
* **Traceback Pass:** Once the peak is found, a localized pass is triggered to map the final alignment string.

## Results
The tool generates **dot plots** that visualize evolutionary distances between species. By dealing with data at the binary level, the aligner provides real biological insights into how plant species have drifted apart over time.

##  Key Takeaways
* **Logic > Language:** Moving from Python to Rust highlighted that software engineering is about problem-solving, not just syntax.
* **Performance:** Bit-packing and memory-efficient algorithms are essential when dealing with 800-million-base genomes.

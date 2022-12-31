use bio::io::bed;
use bio::io::bed::Reader as BReader;
use bio::io::bed::Record as BRecord;
use bio::io::fasta;
use bio::io::fasta::Reader as FReader;
use bio::io::fasta::Record as FRecord;
use bio::pattern_matching::bndm::BNDM;
use mimalloc::MiMalloc;
use nanorand::rand::pcg64::Pcg64;
use nanorand::Rng;
use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::thread;
use std::collections::BinaryHeap;
use std::cmp::Reverse;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;
//i know long struct huh?
struct Names {
    covid: String,
    red: String,
    sed: String,
    hello: String,
    ritchie: String,
    dennis: String,
    bjarne: String,
    idk: String,
    btw: String,
    also: String,
    is: String,
    there1: String,
    more: String,
    words: String,
    ah: String,
    there2: String,
    are: String,
    many: String,
    what: String,
    am: String,
    doing: String,
    with: String,
    my: String,
    life: String,
    oh: String,
    well: String,
    here: String,
    we: String,
    go: String,
    again: String,
    i: String,
    need: String,
    that: String,
    rand: String,
    num: String,
    homie: String,
    duru: String,
    puru: String,
    thread: String,
    read: String,
}
//implementation for Names for genomes
impl Names {
    #[inline (always)]
    fn default() -> Self {
        Names {
            covid: "chr1".to_string(),
            red: "chr2".to_string(),
            sed: "chr3".to_string(),
            hello: "chr4".to_string(),
            ritchie: "chr5".to_string(),
            dennis: "chr6".to_string(),
            bjarne: "chr7".to_string(),
            idk: "chr8".to_string(),
            btw: "chr9".to_string(),
            also: "chr10".to_string(),
            is: "chr11".to_string(),
            there1: "chr12".to_string(),
            more: "chr13".to_string(),
            words: "chr14".to_string(),
            ah: "chr15".to_string(),
            there2: "chr16".to_string(),
            are: "chr17".to_string(),
            many: "chr18".to_string(),
            what: "chr19".to_string(),
            am: "chr20".to_string(),
            doing: "chr22".to_string(),
            with: "chr23".to_string(),
            my: "chr24".to_string(),
            life: "chr25".to_string(),
            oh: "chr26".to_string(),
            well: "chr27".to_string(),
            here: "chr28".to_string(),
            we: "chr29".to_string(),
            go: "chr30".to_string(),
            again: "chr31".to_string(),
            i: "chr32".to_string(),
            need: "chr33".to_string(),
            that: "chr34".to_string(),
            rand: "chr35".to_string(),
            num: "chr36".to_string(),
            homie: "chr37".to_string(),
            duru: "chr38".to_string(),
            puru: "chr39".to_string(),
            thread: "chr40".to_string(),
            read: "chr41".to_string(),
        }
    }
    //puts all strings to a vector, yes im not sane
    #[inline (always)]
    fn into_vector_first(self) -> Vec<String> {
        let a: Vec<String> = vec![
            self.covid,
            self.red,
            self.sed,
            self.hello,
            self.ritchie,
            self.dennis,
            self.bjarne,
            self.idk,
            self.btw,
            self.also,
            self.is,
            self.there1,
            self.more,
            self.words,
            self.ah,
            self.there2,
            self.are,
            self.many,
            self.what,
            self.am,
        ];
        a
    }
    #[inline (always)]
    fn into_vector_second(self) -> Vec<String> {
        let a: Vec<String> = vec![
            self.doing,
            self.with,
            self.my,
            self.life,
            self.oh,
            self.well,
            self.here,
            self.we,
            self.go,
            self.again,
            self.i,
            self.need,
            self.that,
            self.rand,
            self.num,
            self.homie,
            self.duru,
            self.puru,
            self.thread,
            self.read,
        ];
        a
    }
}
#[inline (always)]
pub fn fasta_procedure() {
    thread::spawn(move || {
        //thread 1
        //used to write first 21 entries (counting unwrap_or) to .fa file
        let seed = Names::default();
        let vector = seed.into_vector_first();
        let mut iterate = vector.iter();
        for _ in 0..=20 {
            let mut random = Pcg64::new();
            let mut p_env = env::current_dir().unwrap();
            p_env.push("gen.fa");
            let mut gen_vec: Vec<String> = vec![
                "A".to_string(),
                "C".to_string(),
                "G".to_string(),
                "T".to_string(),
            ];
            gen_vec.extend_from_within(..4);
            gen_vec.extend_from_within(..8);
            gen_vec.extend_from_within(..12);
            gen_vec.extend_from_within(..16);
            random.shuffle(&mut gen_vec);
            let genome = gen_vec.concat().into_bytes();
            let ap_file = File::options().append(true).open(p_env).unwrap();
            let s_buffer = BufWriter::new(ap_file);
            let numbers = random.generate_range(0..=1500);
            let mut id = "id".to_string();
            id.push_str(&numbers.to_string());
            let rec = FRecord::with_attrs(
                iterate.next().unwrap_or(&"chr21".to_string()),
                Some(&id),
                &genome,
            );
            let mut dna_writer = fasta::Writer::new(s_buffer);
            dna_writer.write_record(&rec).unwrap();
            assert_eq!(rec.check(), Ok(()));
        }
    })
    .join()
    .unwrap();
    thread::spawn(|| {
        //thread 2
        //used to write another 21 entries(counting unwrap_or) to .fa file
        let names = Names::default();
        let f_vector = names.into_vector_second();
        let mut it = f_vector.iter();
        for _ in 0..=20 {
            let mut rand = Pcg64::new();
            let mut env = env::current_dir().unwrap();
            env.push("gen.fa");
            let mut gen_vec: Vec<String> = vec![
                "A".to_string(),
                "C".to_string(),
                "G".to_string(),
                "T".to_string(),
            ];
            gen_vec.extend_from_within(..4);
            gen_vec.extend_from_within(..8);
            gen_vec.extend_from_within(..12);
            gen_vec.extend_from_within(..4);
            rand.shuffle(&mut gen_vec);
            let genome = gen_vec.concat().into_bytes();
            let append_file = File::options().append(true).open(env).unwrap();
            let f_buffer = BufWriter::new(append_file);
            let numbers = rand.generate_range(0..=1500);
            let mut id = "id".to_string();
            id.push_str(&numbers.to_string());
            let record = FRecord::with_attrs(
                it.next().unwrap_or(&"chr42".to_string()),
                Some(&id),
                &genome,
            );
            let mut f_writer = fasta::Writer::new(f_buffer);
            f_writer.write_record(&record).unwrap();
            assert_eq!(record.check(), Ok(()));
        }
    })
    .join()
    .unwrap();
    //Reading part and matching
    //2 times recursion
    for _ in 0..=1{
    let mut randomize=Pcg64::new();
    let f_pattern_vec:Vec<&str>=vec!["A","G"];
    let s_pattern_vec:Vec<&str>=vec!["C","T"];
    let mut randomized_vec:Vec<&str>=Vec::with_capacity(3);
    //2 times recursion for pushing
    for _ in 0..=1{
    let pat:usize=randomize.generate_range(0..=1);
    randomized_vec.push(f_pattern_vec.get(pat).unwrap());
    randomized_vec.push(s_pattern_vec.get(pat).unwrap());
    }
    let mut r_env = env::current_dir().unwrap();
    r_env.push("gen.fa");
    let pattern=randomized_vec.concat();
    let read_file = File::open(r_env).unwrap();
    let reader = FReader::new(read_file);
    println!("***FASTA***");
    for record in reader.records() {
        let it = record.unwrap();
        let bndm = BNDM::new(pattern.bytes());
        let occ: Vec<usize> = bndm.find_all(it.seq()).collect();
        if !occ.is_empty() {
            println!("There are {} entries matching {} in {}",occ.len(),pattern,it.desc().unwrap());
        }
    }
}
}
#[inline (always)]
pub fn bed_procedure() {
    //thread 1
    //thread one writes 21 names (with unwrap_or) with randomized positions of genomes with LINE extension and + or - aux
    thread::spawn(move || {
        let names = Names::default();
        let f_vector = names.into_vector_first();
        let mut it = f_vector.iter();
        for _ in 0..=20 {
            //plus minus vector
            let mut pm_vec:Vec<&str>=vec!["+","-"];
            let mut random = Pcg64::new();
            random.shuffle(&mut pm_vec);
            let range:usize=random.generate_range(0..=1);
            let f_rand: u64 = random.generate_range(0..=1000);
            let s_rand: u64 = random.generate_range(1000..=1500);
            let f_score_value:f64=random.generate();
            let s_score_value:i32=random.generate_range(1..=3);
            let actual_score=f64::from(s_score_value)-f_score_value;
            let mut env = env::current_dir().unwrap();
            env.push("dna.bed");
            let f_env = env.clone();
            let s_env = f_env.clone();
            let append_file = File::options().append(true).open(s_env).unwrap();
            let f_buffer = BufWriter::new(append_file);
            let mut record = BRecord::new();
            record.set_chrom(it.next().unwrap_or(&"chr21".to_string()));
            record.set_start(f_rand);
            record.set_end(s_rand);
            record.set_score(&actual_score.to_string());
            record.set_name("LINE");
            record.push_aux(pm_vec.get(range).unwrap());
            let mut f_writer = bed::Writer::new(f_buffer);
            f_writer.write(&record).unwrap();
        }
    })
    .join()
    .unwrap();
    //thread 2
    //does the same as thread 1
    thread::spawn(move || {
        let s_names = Names::default();
        let s_vector = s_names.into_vector_second();
        let mut it = s_vector.iter();
        for _ in 0..=20 {
            let mut pm_vec:Vec<&str>=vec!["+","-"];
            let mut random = Pcg64::new();
            random.shuffle(&mut pm_vec);
            let range:usize=random.generate_range(0..=1);
            let f_rand: u64 = random.generate_range(0..=1000);
            let s_rand: u64 = random.generate_range(1000..=1500);
            let f_score_value:f64=random.generate();
            let s_score_value:i32=random.generate_range(1..=3);
            let actual_score=f64::from(s_score_value)-f_score_value;
            let mut env = env::current_dir().unwrap();
            env.push("dna.bed");
            let f_env = env.clone();
            let s_env = f_env.clone();
            let append_file = File::options().append(true).open(s_env).unwrap();
            let f_buffer = BufWriter::new(append_file);
            let mut record = BRecord::new();
            record.set_chrom(it.next().unwrap_or(&"chr42".to_string()));
            record.set_start(f_rand);
            record.set_end(s_rand);
            record.set_score(&actual_score.to_string());
            record.set_name("LINE");
            record.push_aux(pm_vec.get(range).unwrap());
            let mut f_writer = bed::Writer::new(f_buffer);
            f_writer.write(&record).unwrap();
        }
    })
    .join()
    .unwrap();
    //checks which chrom has the highest length or lowest length
    let mut path = env::current_dir().unwrap();
    path.push("dna.bed");
    let file = File::open(path).unwrap();
    let mut reader = BReader::new(file);
    let mut records_heap:BinaryHeap<(u64,String)>=BinaryHeap::new();
    let mut reverse_heap:BinaryHeap<Reverse<(u64,String)>>=BinaryHeap::new();
    for record in reader.records() {
    let unwrapped=record.unwrap();
    let diff=unwrapped.end()-unwrapped.start();
    let name=unwrapped.chrom().to_string();
    let tup=(diff,name);
    let rev_tup=tup.clone();
    records_heap.push(tup);
    reverse_heap.push(Reverse(rev_tup));
}
//quick doc, i used .0 to remove reverse from tuple
let highest=records_heap.peek().unwrap();
let lowest=reverse_heap.peek().unwrap();
println!("***BED***");
println!("{} has the highest length of {}",highest.1,highest.0);
println!("{} has the lowest length of {}",lowest.0.1,lowest.0.0);
}

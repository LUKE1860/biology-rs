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
            covid: "Covid".to_string(),
            red: "Red".to_string(),
            sed: "Rna".to_string(),
            hello: "Hello".to_string(),
            ritchie: "Ritchie".to_string(),
            dennis: "Dennis".to_string(),
            bjarne: "Bjarne".to_string(),
            idk: "Idk".to_string(),
            btw: "Btw".to_string(),
            also: "Also".to_string(),
            is: "Is".to_string(),
            there1: "There1".to_string(),
            more: "More".to_string(),
            words: "Words".to_string(),
            ah: "Ah".to_string(),
            there2: "There2".to_string(),
            are: "Are".to_string(),
            many: "Many".to_string(),
            what: "What".to_string(),
            am: "Am".to_string(),
            doing: "Doing".to_string(),
            with: "With".to_string(),
            my: "My".to_string(),
            life: "Life".to_string(),
            oh: "Oh".to_string(),
            well: "Well".to_string(),
            here: "Here".to_string(),
            we: "We".to_string(),
            go: "Go".to_string(),
            again: "Again".to_string(),
            i: "I".to_string(),
            need: "Need".to_string(),
            that: "That".to_string(),
            rand: "Rand".to_string(),
            num: "Num".to_string(),
            homie: "Homie".to_string(),
            duru: "Duru".to_string(),
            puru: "Puru".to_string(),
            thread: "Thread".to_string(),
            read: "Read".to_string(),
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
            random.shuffle(&mut gen_vec);
            let genome = gen_vec.concat().into_bytes();
            let ap_file = File::options().append(true).open(p_env).unwrap();
            let s_buffer = BufWriter::new(ap_file);
            let numbers = random.generate_range(0..=1500);
            let mut id = "id".to_string();
            id.push_str(&numbers.to_string());
            let rec = FRecord::with_attrs(
                iterate.next().unwrap_or(&"Position".to_string()),
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
            rand.shuffle(&mut gen_vec);
            let genome = gen_vec.concat().into_bytes();
            let append_file = File::options().append(true).open(env).unwrap();
            let f_buffer = BufWriter::new(append_file);
            let numbers = rand.generate_range(0..=1500);
            let mut id = "id".to_string();
            id.push_str(&numbers.to_string());
            let record = FRecord::with_attrs(
                it.next().unwrap_or(&"Running_low".to_string()),
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
            let mut env = env::current_dir().unwrap();
            env.push("dna.bed");
            let f_env = env.clone();
            let s_env = f_env.clone();
            let append_file = File::options().append(true).open(s_env).unwrap();
            let f_buffer = BufWriter::new(append_file);
            let mut record = BRecord::new();
            record.set_chrom(it.next().unwrap_or(&"Position".to_string()));
            record.set_start(f_rand);
            record.set_end(s_rand);
            record.set_score("0");
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
            let mut env = env::current_dir().unwrap();
            env.push("dna.bed");
            let f_env = env.clone();
            let s_env = f_env.clone();
            let append_file = File::options().append(true).open(s_env).unwrap();
            let f_buffer = BufWriter::new(append_file);
            let mut record = BRecord::new();
            record.set_chrom(it.next().unwrap_or(&"Cordinates".to_string()));
            record.set_start(f_rand);
            record.set_end(s_rand);
            record.set_score("0");
            record.set_name("LINE");
            record.push_aux(pm_vec.get(range).unwrap());
            let mut f_writer = bed::Writer::new(f_buffer);
            f_writer.write(&record).unwrap();
        }
    })
    .join()
    .unwrap();
    //in production
    //implement correction mechanism
    let mut random=Pcg64::new();
    let f_gen:u64=random.generate_range(0..=1500);
    let s_gen:u64=random.generate_range(1500..=1505);
    let mut path = env::current_dir().unwrap();
    path.push("dna.bed");
    let file = File::open(path).unwrap();
    let mut reader = BReader::new(file);
    for record in reader.records() {
    let mut unwrapped=record.unwrap();
    let diff = unwrapped.end() - unwrapped.start();
    if diff <500{
    }
        /*
        let iterated = record.unwrap();
        println!("Length of {} is {}. {:?}", iterated.chrom(), diff,iterated.strand().unwrap())
        */
    }
}

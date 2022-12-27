mod bio;
use std::env;
use std::fs;
use std::fs::File;
use std::path::PathBuf;
//bio::bed_procedure();
//bio::fasta_procedure();
#[inline(always)]
fn create_file(f_path:PathBuf,s_path:PathBuf){
    File::create(f_path).unwrap();
    File::create(s_path).unwrap();
    bio::bed_procedure();
    bio::fasta_procedure();
}
#[inline(always)]
fn truncate_file(f_path:PathBuf,s_path:PathBuf){
    File::options().write(true).truncate(true).open(f_path).unwrap();
    File::options().write(true).truncate(true).open(s_path).unwrap();
}
#[inline(always)]
fn delete_file(f_path:PathBuf,s_path:PathBuf){
fs::remove_file(f_path).unwrap();
fs::remove_file(s_path).unwrap();
}
#[inline(always)]
fn main() {
let mut f_location=env::current_dir().unwrap();
f_location.push("dna.bed");
let mut s_location=env::current_dir().unwrap();
s_location.push("gen.fa");
let empty=String::from("");
let args:Vec<String>=env::args().collect();
let arg=args.get(1).unwrap_or(&empty);
match arg.trim(){
"-c"=>create_file(f_location,s_location),
"-t"=>truncate_file(f_location,s_location),
"-d"=>delete_file(f_location,s_location),
"-h"=>println!("-c creates files and runs a program.\n 
-t truncates files and doesn't run a program.\n 
-d deletes all files and doesn't run a program (i wouldn't let you).\n
without arguments program runs as normal.\n
If you are running the program first time use -c."),
_=>{
bio::bed_procedure();
bio::fasta_procedure();
},
};
}

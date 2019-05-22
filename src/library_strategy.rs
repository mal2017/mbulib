#[derive(Debug)]
pub enum LibraryType {
    R1Sense,
    R2Sense,
    Unstranded,
    ATAC,
    DNASE,

}

#[derive(Debug)]
pub enum CountStrand {
    Plus,
    Minus,
    Both,
}

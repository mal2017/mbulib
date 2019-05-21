quick_error! {
    #[derive(Debug, Clone)]
    pub enum InvalidRecordError {
        UnmappedRecord {
            description("Record is unmapped")
        }
        UnknownScaffold {
            description("Record tid not present in provided scaffold dict")
        }
        UnpairedRecord {
            description("Record is unpaired but your library interpretation requires pairs")
        }
    }
}

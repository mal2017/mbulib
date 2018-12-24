use rust_htslib::bam;
use std::collections::HashMap;
use linear_map::LinearMap;
use std::string::String;

// TODO handle panics
/// From a HeaderView return a new Header with a the SO tag changed.
pub fn edit_hdr_srt_tag(hdrv: &bam::HeaderView, new_srt_tag: &str) -> bam::header::Header {
	let hdr = bam::Header::from_template(hdrv);
	let mut hm = hdr.to_hashmap().to_owned();
	let hd = &hm.get("HD").unwrap().to_owned();
	// TODO Assert only 1 HD entry
	let mut record = hd[0].to_owned();
	&mut record.insert("SO".to_string(),new_srt_tag.to_string());
	hm.insert("HD".to_owned(),vec![record]);
	from_hashmap(&hm)
}

/// Get the value of the Header
pub fn get_hdr_srt_tag(hdrv: &bam::HeaderView) -> Option<String> {
    let hdr = bam::Header::from_template(hdrv);
    let mut hm = hdr.to_hashmap().to_owned();
    let hd = &hm.get("HD").unwrap().to_owned();
    let mut record = hd[0].to_owned();
    match record.get(&"SO".to_string()) {
        Some(s) => Some(s.to_owned()),
        None => None,
    }
}


// TODO Handle panics
/// From a bam header stored as a HashMap storing a vector of LinearMaps, spit out a new header.
fn from_hashmap(hm: &HashMap<String, Vec<LinearMap<String, String>>>) ->  bam::header::Header {
	let tags = vec!["HD","SQ","PG","CL","RG"];
	let bam_header_strings = tags.iter()
								 .map(|a| hdr_tag_to_string_vector(a, &hm));
	let mut hdr_records_string: Vec<String> = vec![];
	for x in bam_header_strings {
		match x {
			Some(x) => {
				let mut y = x.join("\n");
				hdr_records_string.push(y);
			},
			None => (),
		}
	};
	let new_hdr_view = bam::HeaderView::from_bytes(hdr_records_string.join("\n").as_bytes());
	bam::Header::from_template(&new_hdr_view)
}

// TODO hanlde panics
/// From a Linear Map generates a string in the format "KEY:VALUE".
/// Used for working with bam headers stored in HashMaps storing vectors of LinearMaps,
/// as is currently an option in rust_htslib.
fn linear_map_to_hdr_entry_string(lm: &LinearMap<String, String>, tag: &str) -> String {
	let entries: Vec<String> = lm.iter()
								 .map(|a| format!("{}:{}",a.0,a.1))
								 .collect();
	format!("@{}\t{}",tag,entries.join("\t"))
}

/// For a supplied header tag and a ref to a hashmap of header data,
/// returns a vector of string representations.
fn hdr_tag_to_string_vector(tag: &str, hm: &HashMap<String, Vec<LinearMap<String, String>>>) -> Option<Vec<String>> {
	let record = hm.get(tag);
	match record {
		Some(r) => {
			let string_records: Vec<String> = r.iter()
		  									   .map(|a| linear_map_to_hdr_entry_string(a, tag))
		  									   .collect();
		  	Some(string_records)
		},
		None => None,
	}
}

#[cfg(test)]
mod tests {
    use crate::bam::header::*;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use rust_htslib::bam::*;
    use std::path::Path;

    #[test]
    fn grab_header_works() {
        let bampath = Path::new("test/aln.mini.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let so_tag = get_hdr_srt_tag(bam.header()).unwrap();
        assert_eq!(so_tag,"coordinate");
    }

    #[test]
    fn edit_header_works() {
        let bampath = Path::new("test/aln.mini.bam");
        let mut bam = bam::Reader::from_path(bampath).unwrap();
        let new_hdrv = HeaderView::from_header(&edit_hdr_srt_tag(bam.header(), "test"));
        let so_tag = get_hdr_srt_tag(&new_hdrv).unwrap();
        assert_eq!(so_tag,"test");
    }
}

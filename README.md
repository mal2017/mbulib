# rqmap

## Background

The main goal of this crate is to provide ways to easily query the mapped
locations of large sets of NGS reads.

This is a pretty standard task for many NGS analyses, but
I found myself writing a lot of boilerplate to do this in some of my personal
crates and cli tools.

## Major goals:

* easily ingest a bam and create an RQMap (wrapper around bio::AnnotMap)
* easily filter or preprocess the reads during this process
* easily query the ingested data by overlaps with bio_types::Contig


## TODOs

 - read level filter option for map construction
 - Specify fragments or reads to be added to the map
 - Single high level function/macro for creating an RQMap
 - Single high level function/macro for querying the overlaps

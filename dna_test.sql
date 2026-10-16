DROP TABLE d;
DROP EXTENSION dna;
CREATE EXTENSION dna;
CREATE TABLE d (id SERIAL PRIMARY KEY, dn DNA);
INSERT INTO d(dn) VALUES
 ('GCGTC'), ('ACGCGTC'),('aacgtagca'),('a'),('gggctcgag'),('ATGC'), ('GCGTC'), ('ACGCGTC'), ('ACGTAATTT'), ('AAGTACGCGCG'), ('ATGTACC'),('AA');
SELECT * FROM d WHERE length(dn) = 3 ;

DROP TABLE d;
DROP EXTENSION dna;
CREATE EXTENSION dna;
CREATE TABLE d (id SERIAL PRIMARY KEY, dn KMER);

INSERT INTO d(dn) VALUES
 ('ATGC'), ('GCGTC'), ('ACGCGTC'), ('ACGTAATTT'), ('AAGTACGCGCG'), ('ATGTACC'),('AA');

SELECT * FROM d WHERE contains('ANGTA', dn);
SELECT * FROM d WHERE starts_with('ACG', dn);
SELECT k.kmer FROM generate_kmers('ACGTACGT', 6) AS k(kmer);
SELECT * FROM d WHERE dn = 'GCGTC';
SELECT * FROM d WHERE length (dn)= 5;

-- Test less than
SELECT * FROM d WHERE dn < 'GCTA'::kmer;
-- Test greater than
SELECT * FROM d WHERE dn > 'ACT'::kmer;
-- Test range queries
SELECT * FROM d WHERE dn BETWEEN 'ACT'::kmer AND 'ACTGA'::kmer;
-- Test sorting
SELECT * FROM d ORDER BY dn;

SELECT COUNT(*) FROM d;

CREATE INDEX d_kmer_index ON d USING spgist (dn);
SET enable_seqscan = off;
EXPLAIN ANALYZE SELECT * FROM d WHERE dn = 'ATGC';
EXPLAIN ANALYZE SELECT * FROM d WHERE 'AT' @> dn;

DROP TABLE d;
DROP EXTENSION dna;
#!/usr/bin/env python3

import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FIXTURES = Path(os.environ.get("NFSCRNA_FIXTURES", "/tmp/nfscrna-gex-fixtures"))
R_IMAGE = "ghcr.io/obenauflab/catchseq:latest"
NEXTFLOW = os.environ.get("NEXTFLOW_BIN") or shutil.which("nextflow")
LIBRARIES = ("S1_ctrl_1", "S2_ctrl_2", "S3_trt_1", "S4_trt_2")


def docker_image_available(image):
    try:
        result = subprocess.run(
            ["docker", "image", "inspect", image], capture_output=True, text=True
        )
    except FileNotFoundError:
        return False
    return result.returncode == 0


def generate_fixtures(image):
    FIXTURES.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-u",
            f"{os.getuid()}:{os.getgid()}",
            "-e",
            f"NFSCRNA_FIXTURES=/fixtures/{FIXTURES.name}",
            "-e",
            "HOME=/tmp",
            "-v",
            f"{REPO}:/repo:ro",
            "-v",
            f"{FIXTURES.parent}:/fixtures",
            image,
            "Rscript",
            "--vanilla",
            "/repo/tests/make_fixtures.R",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"fixture generation failed: {result.stdout[-4000:]}{result.stderr[-4000:]}"
        )
    return FIXTURES


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tree_hashes(root):
    return {
        str(path.relative_to(root)): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def h5_columns(path, group):
    result = subprocess.run(
        ["h5ls", f"{path}/{group}"], capture_output=True, text=True
    )
    if result.returncode != 0:
        raise AssertionError(f"h5ls failed on {path}/{group}: {result.stderr}")
    return [line.split()[0] for line in result.stdout.splitlines() if line.strip()]


def pdf_pages(path):
    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(f"pdfinfo failed on {path}: {result.stderr}")
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split()[1])
    raise AssertionError(f"no page count in pdfinfo output for {path}")


class WorkflowSmokeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not NEXTFLOW:
            raise unittest.SkipTest("nextflow is not available")
        if not docker_image_available(R_IMAGE):
            raise unittest.SkipTest(f"{R_IMAGE} is not available")
        if not FIXTURES.is_dir():
            generate_fixtures(R_IMAGE)
        cls.root = Path(tempfile.mkdtemp(prefix="nfscrna-workflow-", dir="/tmp"))
        cls.launch = cls.root / "launch"
        cls.launch.mkdir()
        cls.work = cls.root / "work"
        cls.precomputed = cls.root / "precomputed"
        for library in LIBRARIES:
            for matrix in ("raw_feature_bc_matrix", "filtered_feature_bc_matrix"):
                target = cls.precomputed / library / matrix
                target.mkdir(parents=True)
                for name in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
                    shutil.copyfile(FIXTURES / library / name, target / name)
        cls.annotation = cls.root / "annotation.gtf"
        shutil.copyfile(FIXTURES / "annotation.gtf", cls.annotation)
        cls.config = cls.launch / "smoke.config"
        cls.config.write_text(
            "docker {\n"
            "    enabled = true\n"
            "    remove = true\n"
            "    runOptions = '-u $(id -u):$(id -g)'\n"
            "}\n"
            "\n"
            "process {\n"
            "    cpus = 1\n"
            "    memory = 6.GB\n"
            "    time = 10.min\n"
            "    maxForks = 1\n"
            "    withName: preprocessSingleCellData {\n"
            "        cpus = 1\n"
            "        memory = 6.GB\n"
            "        time = 10.min\n"
            "    }\n"
            "    withName: createOverviewPlots {\n"
            "        cpus = 1\n"
            "        memory = 6.GB\n"
            "        time = 10.min\n"
            "    }\n"
            "    withName: identifyDEGenes {\n"
            "        cpus = 1\n"
            "        memory = 6.GB\n"
            "        time = 10.min\n"
            "    }\n"
            "}\n"
            "\n"
            "executor {\n"
            "    $local {\n"
            "        cpus = 1\n"
            "        memory = '6 GB'\n"
            "    }\n"
            "}\n"
            "\n"
            "env {\n"
            "    OPENBLAS_NUM_THREADS = '1'\n"
            "    OMP_NUM_THREADS = '1'\n"
            "    MKL_NUM_THREADS = '1'\n"
            "    HOME = '/tmp'\n"
            "}\n"
        )
        cls.before = tree_hashes(FIXTURES)
        cls.full = cls.run_workflow(
            "full",
            cls.samplesheet(cls.root / "full.csv", LIBRARIES),
            [
                "--filter",
                "false",
                "--output_seurat",
                "true",
                "--output_sce",
                "true",
                "--output_anndata",
                "true",
                "--run_sct",
                "true",
                "--min_detected_features",
                "5",
                "--max_mt_percent",
                "100",
                "--minReads",
                "1",
            ],
        )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)
        if cls.root.exists():
            subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "-v",
                    f"{cls.root.parent}:/cleanup",
                    "busybox",
                    "rm",
                    "-rf",
                    f"/cleanup/{cls.root.name}",
                ],
                capture_output=True,
            )
            shutil.rmtree(cls.root, ignore_errors=True)

    @classmethod
    def samplesheet(cls, path, libraries):
        rows = [
            "SampleName,Condition,Replicate,Lane,LibraryType,R1,R2,CellNumber,Chemistry"
        ]
        for library in libraries:
            sample, condition, replicate = library.split("_")
            rows.append(
                f"{sample},{condition},{replicate},1,CellRanger,"
                f"{cls.precomputed / library},,NA,10X"
            )
        path.write_text("\n".join(rows) + "\n")
        return path

    @classmethod
    def run_workflow(cls, name, samplesheet, extra, libraries=True, de_rds=None):
        output = cls.launch / f"out_{name}"
        command = [
            NEXTFLOW,
            "-log",
            str(cls.launch / f"nextflow_{name}.log"),
            "-c",
            str(cls.config),
            "run",
            str(REPO / "nextflow" / "main.nf"),
            "--outputDir",
            str(output),
            "--reportsDir",
            str(cls.launch / f"reports_{name}"),
            "--mapper",
            "CellRanger",
            "--annotation",
            str(cls.annotation),
            "--baseline",
            "ctrl",
            "-work-dir",
            str(cls.work),
            "-ansi-log",
            "false",
        ]
        if libraries:
            command.extend(["--libraries", str(samplesheet)])
        if de_rds is not None:
            command.extend(["--de_rds", str(de_rds)])
        command.extend(str(argument) for argument in extra)
        result = subprocess.run(
            command, capture_output=True, text=True, cwd=cls.launch
        )
        return result, output

    def test_full_workflow_outputs(self):
        result, output = self.full
        self.assertEqual(result.returncode, 0, result.stdout[-6000:] + result.stderr[-6000:])

        for library in LIBRARIES:
            tag = library
            for name in (
                f"SCE/filtered/{tag}_nf-scRNASeq_filtered_seurat_sce.rds.gz",
                f"SCE/filtered/{tag}_nf-scRNASeq_filtered_sce.rds.gz",
                f"SCE/raw/{tag}_nf-scRNASeq_unfiltered_seurat_sce.rds.gz",
                f"SCE/raw/{tag}_nf-scRNASeq_unfiltered_sce.rds.gz",
                f"AnnData/filtered/{tag}_nf-scRNASeq_filtered.h5ad",
                f"AnnData/raw/{tag}_nf-scRNASeq_unfiltered.h5ad",
                f"Plots/Overview/{tag}/{tag}_nf-scRNASeq_overview.pdf",
                f"Plots/Overview/{tag}/{tag}_nf-scRNASeq_overview_QC_summary.tsv",
                f"Plots/Overview/{tag}/{tag}_nf-scRNASeq_overview_cell_metadata.tsv",
            ):
                path = output / name
                self.assertTrue(path.is_file(), f"missing {name}: {sorted(p.name for p in (output).rglob('*'))[:40]}")
                self.assertGreater(path.stat().st_size, 0, name)

        overview = output / "Plots/Overview/S1_ctrl_1/S1_ctrl_1_nf-scRNASeq_overview.pdf"
        self.assertGreater(pdf_pages(overview), 0)
        qc_pdf = output / "Plots/Overview/S1_ctrl_1/S1_ctrl_1_nf-scRNASeq_overview_QC_MT_content_filtered.pdf"
        self.assertGreater(pdf_pages(qc_pdf), 0)

        summary = (
            output / "Plots/Overview/S1_ctrl_1/S1_ctrl_1_nf-scRNASeq_overview_QC_summary.tsv"
        ).read_text().splitlines()
        self.assertEqual(len(summary), 2, summary)
        self.assertEqual(summary[1].split("\t")[1], "ctrl")

        for library in LIBRARIES:
            sample, condition, replicate = library.split("_")
            lines = (
                output
                / f"Plots/Overview/{library}/{library}_nf-scRNASeq_overview_QC_summary.tsv"
            ).read_text().splitlines()
            self.assertEqual(len(lines), 2, lines)
            columns = lines[0].split("\t")
            values = lines[1].split("\t")
            self.assertEqual(values[columns.index("Sample")], sample)
            self.assertEqual(values[columns.index("Condition")], condition)
            self.assertEqual(values[columns.index("Replicate")], f"{condition}_{replicate}")
            self.assertGreater(int(values[columns.index("n_cells")]), 0)

        de_all = output / "DE/GENES/DE_DESEQ2_trt-vs-ctrl_GENES_all.tsv.gz"
        self.assertTrue(de_all.is_file(), sorted(p.name for p in (output / "DE").rglob("*")))
        with gzip.open(de_all, "rt") as handle:
            self.assertGreater(len(handle.read().splitlines()), 1)
        self.assertTrue((output / "DE/GENES/allSamples_nf-scRNASeq_DE_summary.rds.gz").is_file())
        skipped = (output / "DE/GENES/allSamples_nf-scRNASeq_DE_skipped.tsv").read_text().splitlines()
        self.assertEqual(len(skipped), 1, skipped)
        self.assertTrue(list((output / "DE/Plots").glob("*.pdf")))

        obs = h5_columns(
            output / "AnnData/filtered/S1_ctrl_1_nf-scRNASeq_filtered.h5ad", "obs"
        )
        self.assertIn("Sample", obs)
        self.assertIn("Condition", obs)
        self.assertIn("Replicate", obs)
        self.assertFalse([column for column in obs if "CaTCH" in column], obs)
        self.assertEqual(tree_hashes(FIXTURES), self.before)

    def test_flags_off_keep_internal_seurat(self):
        samplesheet = self.samplesheet(self.root / "single.csv", LIBRARIES[:1])
        result, output = self.run_workflow(
            "flags_off",
            samplesheet,
            [
                "--filter",
                "false",
                "--output_seurat",
                "false",
                "--output_sce",
                "false",
                "--output_anndata",
                "false",
                "--run_sct",
                "false",
                "--min_detected_features",
                "5",
                "--max_mt_percent",
                "100",
                "--minReads",
                "1",
            ],
        )
        self.assertEqual(result.returncode, 0, result.stdout[-6000:] + result.stderr[-6000:])
        self.assertFalse((output / "SCE").exists(), sorted(p.name for p in output.rglob("*")))
        self.assertFalse((output / "AnnData").exists())
        internal = list(self.work.rglob("*_filtered_seurat_sce.rds.gz"))
        self.assertTrue(internal, "filtered Seurat object was not produced internally")
        self.assertTrue(list((output / "Plots/Overview").rglob("*_overview.pdf")))
        self.assertTrue((output / "DE/GENES").is_dir())
        skipped = (
            output / "DE/GENES/allSamples_nf-scRNASeq_DE_skipped.tsv"
        ).read_text().splitlines()
        self.assertEqual(len(skipped), 2, skipped)
        self.assertIn("single_condition", skipped[1])

    def test_de_rds_only(self):
        rds = (
            self.full[1]
            / "SCE/filtered/S1_ctrl_1_nf-scRNASeq_filtered_seurat_sce.rds.gz"
        )
        self.assertTrue(rds.is_file(), rds)
        result, output = self.run_workflow(
            "de_rds",
            None,
            ["--filter", "false"],
            libraries=False,
            de_rds=rds,
        )
        self.assertEqual(result.returncode, 0, result.stdout[-6000:] + result.stderr[-6000:])
        self.assertTrue((output / "DE/GENES/de_rds_nf-scRNASeq_DE_summary.rds.gz").is_file(),
            sorted(p.name for p in (output / "DE").rglob("*")))
        self.assertFalse((output / "SCE").exists())
        self.assertFalse((output / "RAWMAPPING").exists())
        log = (self.launch / "nextflow_de_rds.log").read_text()
        submitted = [line for line in log.splitlines() if "Submitted process" in line]
        self.assertTrue(submitted, log[-4000:])
        for process in ("useCellrangerData", "preprocessSingleCellData"):
            self.assertFalse(
                [line for line in submitted if process in line],
                f"{process} was submitted in the de_rds run: {submitted}",
            )


if __name__ == "__main__":
    unittest.main()

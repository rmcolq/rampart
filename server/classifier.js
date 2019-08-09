const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');
const Deque = require("collections/deque");
const { getReadTime } = require('./readTimes');
const { prettyPath, warn, verbose } = require('./utils');

/**
 * This file defines the global.classifyingQueue handler, which processes demuxed
 * FASTQ files. Note that the data has already been partially extracted
 * from these files, and a corresponding entry is present in global.datastore
 *
 * global.classifyingQueue contains objets of `[pointer, fileToClassify]` where
 * the pointer defines the idx of the partial info in global.datastore.
 * global.classifyingQueue is a `Deque`
 *
 * The interface between node and the classifier (currently a snakemake which
 * calls kraken2) is contained here
 */
const classifyingQueue = new Deque();

classifyingQueue.observeRangeChange(() => {classifier();});
const addToClassifyingQueue = (thing) => classifyingQueue.push(thing);

const call_python_bin_finder = (fastq) => new Promise((resolve, reject) => {
    const basename = path.basename(fastq, ".fastq");
    const workDir = path.dirname(global.config.demuxedPath);
    const snakefile = path.join(__dirname, "..", "snakefiles/kraken/Snakefile");
    const pyprog = spawn('python3', [
        "./binlorry/binlorry/get_key_values.py",
        fastq,
        "barcode"
    ]);
    let stdout = "";
    let stderr = "";
    pyprog.stdout.on('data', (data) => {stdout+=data});
    pyprog.stderr.on('data', (data) => {stderr+=data});

    // stochastically mock failure
    if (global.MOCK_FAILURES && Math.random() < 0.05) {
        reject("Mock classifying failure")
    }

    pyprog.on('close', (code) => {
        if (code === 0) {
            if (stderr) console.log(stderr);
            resolve(stdout);
        } else {
            reject(stderr)
        }
    });
});

const call_snakemake_classifier = (fastq, barcodes) => new Promise((resolve, reject) => {
    const basename = path.basename(fastq, ".fastq");
    const workDir = path.dirname(global.config.demuxedPath);
    const snakefile = path.join(__dirname, "..", "snakefiles/kraken/Snakefile");
    const pyprog = spawn('snakemake', [
        "--snakefile", snakefile,
        "--config",
        "workdir=" + workDir,
        "sample=" + basename,
        "barcodes=" + barcodes
    ]);
    let stdout = "";
    let stderr = "";
    pyprog.stdout.on('data', (data) => {stdout+=data});
    pyprog.stderr.on('data', (data) => {stderr+=data});

    // stochastically mock failure
    if (global.MOCK_FAILURES && Math.random() < 0.05) {
        reject("Mock classifying failure")
    }

    pyprog.on('close', (code) => {
        if (code === 0) {
            if (stderr) console.log(stderr);
            resolve(code);
        } else {
            reject(stderr)
        }
    });
});

let isRunning = false; // only want one classifying thread at a time!
const classifier = async () => {

    if (isRunning) {
        verbose("[classifier] called but already running");
        return;
    }

    if (classifyingQueue.length) {
        isRunning = true;
        let results;
        const [datastoreAddress, fileToClassify] = classifyingQueue.shift();
        try {
            verbose(`[classifier] queue length: ${classifyingQueue.length+1}. classifying ${prettyPath(fileToClassify)}`);
            barcodes = await call_python_bin_finder(fileToClassify);
            results = await call_snakemake_classifier(fileToClassify, barcodes);
            //global.datastore.addClassifiedFastq(datastoreAddress, results);
            verbose(`[classifier] Classified ${prettyPath(fileToClassify)}. Read time: ${getReadTime(fileToClassify)}.`);
        } catch (err) {
            console.trace(err);
            //warn(`classifying ${fileToClassify.split("/").slice(-1)[0]}: ${err}`);
        }
        isRunning = false;
        classifier(); // recurse
    }
}

module.exports = {
    classifier,
    addToClassifyingQueue,
}

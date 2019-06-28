const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');
const Deque = require("collections/deque");
const { getReadTime } = require('./readTimes');
const { prettyPath, warn, verbose } = require('./utils');

/**
 * This file defines the global.refFinderQueue handler, which processes demuxed
 * FASTQ files. Note that the data has already been partially extracted
 * from these files, and a corresponding entry is present in global.datastore
 *
 * global.refFinderQueue contains objects of `[pointer, fileToClassify]` where
 * the pointer defines the idx of the partial info in global.datastore.
 * global.refFinderQueue is a `Deque`
 *
 * The interface between node and the ref finding process (currently a python script which
 * calls kraken2) is contained here
 */
const refFinderQueue = new Deque();

refFinderQueue.observeRangeChange(() => {ref_finder();});
const addToRefFinderQueue = (thing) => refFinderQueue.push(thing);

/*const call_python_ref_finder = (fastq) => new Promise((resolve, reject) => {
    const pyprog = spawn('python3', [
        "./server/map_single_fastq.py",
        "-c", global.config.coordinateReferencePath,
        "-p", global.config.referencePanelPath,
        "-f", fastq
    ]);
    let stdout = "";
    let stderr = "";
    pyprog.stdout.on('data', (data) => {stdout+=data});
    pyprog.stderr.on('data', (data) => {stderr+=data});

    // stochastically mock failure
    if (global.MOCK_FAILURES && Math.random() < 0.05) {
        reject("Mock mapping failure")
    }

    pyprog.on('close', (code) => {
        // console.log(`Python script finished. Exit code ${code}`);
        if (code === 0) {
            if (stderr) console.log(stderr);
            resolve(JSON.parse(stdout));
        } else {
            reject(stderr)
        }
    });
});*/

let isRunning = false; // only want one mapping thread at a time!
const ref_finder = async () => {

    /* the ref_finder can _only_ run _if_ we have defined both a reference panel and a
    "main" reference config. (Perhaps this could be relaxed in the future) */
    /*if (!(global.config.reference && global.config.referencePanel.length)) {
        verbose(`Cannot start ref_finder without main reference (provided: ${!!global.config.reference}) AND reference panel (provided: ${!!global.config.referencePanel.length})`);
        return;
    }*/

    if (isRunning) {
        verbose("[ref_finder] called but already running");
        return;
    }

    if (refFinderQueue.length) {
        isRunning = true;
        let results;
        const [datastoreAddress, fileToClassify] = refFinderQueue.shift();
        try {
            verbose(`[ref_finder] queue length: ${refFinderQueue.length+1}. Finding reference for ${prettyPath(fileToClassify)}`);
            //results = await call_python_ref_finder(fileToClassify);
            //global.datastore.addMappedFastq(datastoreAddress, results);
            verbose(`[ref_finder] Found reference for ${prettyPath(fileToClassify)}. Read time: ${getReadTime(fileToClassify)}.`);
        } catch (err) {
            console.trace(err);
            warn(`Finding reference for ${fileToClassify.split("/").slice(-1)[0]}: ${err}`);
        }
        isRunning = false;
        ref_finder(); // recurse
    }
}

module.exports = {ref_finder, addToRefFinderQueue};

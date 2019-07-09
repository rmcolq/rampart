const fs = require('fs')
const path = require('path')
const { promisify } = require('util');
const { addToRefFinderQueue } = require("./refFinder");
const { save_coordinate_reference_as_fasta, addToMappingQueue } = require("./mapper");
const { addToDemuxQueue } = require("./demuxer");
const readdir = promisify(fs.readdir);
const { setReadTime, getReadTime, setEpochOffset } = require('./readTimes');
const { prettyPath, verbose, log, deleteFolderRecursive } = require('./utils');

const getFastqsFromDirectory = async (dir) => {
    let fastqs = (await readdir(dir))
        .filter((j) => j.endsWith(".fastq"))
        .sort((a, b) => parseInt(a.match(/\d+/), 10) > parseInt(b.match(/\d+/), 10) ? 1 : -1)
        .map((j) => path.join(dir, j));
    return fastqs;
}

const startUp = async ({emptyDemuxed=false}={}) => {

  log("Rampart Start Up")
  if (!fs.existsSync(global.config.basecalledPath) || !fs.existsSync(global.config.demuxedPath)) {
    verbose("[startUp] no basecalled dir / demuxed dir");
    return false;
  }

  /* Create & clear a temporary directory to store data in during the run */
  if (fs.existsSync(global.config.rampartTmpDir)) {
    deleteFolderRecursive(global.config.rampartTmpDir);
  }
  fs.mkdirSync(global.config.rampartTmpDir);

  if (global.config.reference) {
    /* the python mapping script needs a FASTA of the main reference */
    global.config.coordinateReferencePath = save_coordinate_reference_as_fasta(global.config.reference.sequence, global.config.rampartTmpDir);
  }

  /* LIST BASECALLED & DEMUXED FILES */
  const unsortedBasecalledFastqs = await getFastqsFromDirectory(global.config.basecalledPath);
  log(`  * Found ${unsortedBasecalledFastqs.length} basecalled fastqs in ${prettyPath(global.config.basecalledPath)}`);

  let unsortedDemuxedFastqs = [];
  if (emptyDemuxed) {
    log(`  * Clearing the demuxed folder (${prettyPath(global.config.demuxedPath)})`);
    const demuxedFilesToDelete = await readdir(global.config.demuxedPath);
    for (const file of demuxedFilesToDelete) {
      fs.unlinkSync(path.join(global.config.demuxedPath, file));
    }
  } else {
    unsortedDemuxedFastqs = await getFastqsFromDirectory(global.config.demuxedPath);
    log(`  * Found ${unsortedDemuxedFastqs.length} basecalled fastqs in ${prettyPath(global.config.demuxedPath)}`);
  }

  /* EXTRACT TIMESTAMPS FROM THOSE FASTQs */
  log(`  * Getting read times from basecalled and demuxed files`);
  for (let fastq of unsortedBasecalledFastqs) {
    await setReadTime(fastq);
  }
  for (let fastq of unsortedDemuxedFastqs) {
    await setReadTime(fastq);
  }
  setEpochOffset();



  /* sort the fastqs via these timestamps and push onto the appropriate deques */
  log(`  * Sorting available basecalled and demuxed files based on read times`);
  unsortedDemuxedFastqs
    .sort((a, b) => getReadTime(a)>getReadTime(b) ? 1 : -1)
    .forEach((f) => {
      addToMappingQueue(f);
      addToRefFinderQueue(f);
      global.fastqsSeen.add(path.basename(f));
    })

  const demuxedBasenames = unsortedDemuxedFastqs.map((name) => path.basename(name))
  unsortedBasecalledFastqs
    .filter((fastqPath) => !demuxedBasenames.includes(path.basename(fastqPath)))
    .sort((a, b) => getReadTime(a)>getReadTime(b) ? 1 : -1)
    .forEach((f) => {
      addToDemuxQueue(f);
      global.fastqsSeen.add(path.basename(f))
    });

  log`RAMPART start up FINISHED\n`;
  return true;
}

module.exports = {startUp}

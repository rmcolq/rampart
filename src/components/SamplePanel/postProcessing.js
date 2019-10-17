import React, {useState} from 'react';
import PropTypes from "prop-types";
import Modal from "../modal";
import { IoMdPlay } from "react-icons/io";

const ChooseMinReadLength = ({userSettings, setUserSettings}) => (
    <>
        <h4>Minimum read length:</h4>
        <input type="text" value={userSettings.min_length} onChange={(e) => setUserSettings({...userSettings, min_length: e.target.value})}/>
    </>
);

const ChooseMaxReadLength = ({userSettings, setUserSettings}) => (
    <>
        <h4>Maximum read length:</h4>
        <input type="text" value={userSettings.max_length} onChange={(e) => setUserSettings({...userSettings, max_length: e.target.value})}/>
    </>
);


export const getPostProcessingMenuItems = (config, setPostProcessingState) => {
    return config.pipelines.processing
        .filter((pipeline) => pipeline.run_per_sample)
        .map((obj) => ({
            label: obj.name,
            callback: () => setPostProcessingState(obj)
        }));
};

const createInitialState = (pipeline) => {
    console.log("creating initial state");
    const initialState = {};
    if (pipeline.options.min_length) initialState.min_length = 0; // TODO -- get min read length from dataset
    if (pipeline.options.max_length) initialState.max_length = 100000; // TODO -- get max read length from dataset
    return initialState;
}

export const PostProcessingRunner = ({pipeline, dismissModal, socket, sampleName}) => {
    const [userSettings, setUserSettings] = useState(() => createInitialState(pipeline));

    const send = () => {
        console.log("triggerPostProcessing")
        socket.emit('triggerPostProcessing', {pipeline, sampleName, userSettings});
        dismissModal();
    };

    return (
        <Modal dismissModal={dismissModal}>
            <h2>{pipeline.name}</h2>
            {Object.keys(userSettings).length ? (
                <>
                    <h3>This pipeline requests the following options:</h3>
                    {(userSettings.min_length !== undefined) ? (
                        <ChooseMinReadLength userSettings={userSettings} setUserSettings={setUserSettings}/>
                    ) : null}
                    {(userSettings.max_length !== undefined) ? (
                        <ChooseMaxReadLength userSettings={userSettings} setUserSettings={setUserSettings}/>
                    ) : null}
                </>
            ) : null }
            <button className="modernButton" onClick={send}>
                <div><IoMdPlay/><span>TRIGGER</span></div>
            </button>
        </Modal>
    )
};

PostProcessingRunner.propTypes = {
    dismissModal: PropTypes.func.isRequired,
    sampleName: PropTypes.string.isRequired,
    pipeline: PropTypes.oneOfType([PropTypes.bool, PropTypes.object]).isRequired,
    socket: PropTypes.object.isRequired
};
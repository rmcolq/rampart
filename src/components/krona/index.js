import React, {useState} from 'react';
import InfoRowKrona from "./infoRowKrona";
import {ContextMenuTrigger} from "react-contextmenu";
import {IconContext} from "react-icons";
import {MdReorder} from "react-icons/md";
import Menu from "../SamplePanel/menu";

const Krona = ({config, viewOptions, canExpand, socket}) => {

    /* -----------    STATE MANAGEMENT    ------------------- */
    const [expanded, setExpanded] = useState(false);
    const [singleRow, setSingleRow] = useState(true);
    const [showSinglePanel, setShowSinglePanel] = useState(true);
    const [transitionInProgress, setTransitionInProgress] = useState(false);
    const transitionStarted = (duration = 600) => { /* CSS transition is 0.5s */
        setTransitionInProgress(true);
        setTimeout(() => setTransitionInProgress(false), duration);
    }
    const toggleExpanded = () => {
        transitionStarted();
        setExpanded(!expanded);
    }

    /* ------------- MENU OPTIONS -------------------- */
    const menuItems = [];
    if (!expanded) {
        menuItems.push({label: "Expand panel", callback: toggleExpanded})
    } else {
        menuItems.push({label: "Contract panel", callback: toggleExpanded})
        menuItems.push({label: "Show All (horisontally)", callback: () => {transitionStarted(); setShowSinglePanel(false); setSingleRow(true);}});
    }


    function load_home() {
        document.getElementById("content").innerHTML = '<object type="text/html" data="krona.html" ></object>';
    }

    /* ---------------   WHAT CHARTS DO WE RENDER?   -------------- */
    /*const renderCharts = () => {
        if (!expanded) return null;
        return <div className="panelFlexColumn">
            <div className="panelFlexRow">
                {load_home()}
            </div>
        </div>;
    }




    {transitionInProgress ? null : renderCharts()}*/

    /* ----------------- R E N D E R ---------------- */
    const sampleColour = 'rgb(220,220,220)';
    return (
        <div
            className={`panelContainer ${expanded ? "expanded" : "collapsed"} ${singleRow ? "singleRow" : "multiRow"}`}
            style={{borderColor: sampleColour}}
        >
            <InfoRowKrona
                sampleColour={sampleColour}
                enableUserInteraction={canExpand}
                menuItems={menuItems}
                handleClick={toggleExpanded}
                isExpanded={expanded}
            />
        </div>
    );
};

export default Krona;


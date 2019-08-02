import React, {useState} from 'react';
import InfoRowKrona from "./infoRowKrona";

/**
 * Why are we using this transition / setTimeout stuff?
 *    The charts, upon initial rendering, calculate the SVG dimentions from the DOM they're in.
 *    Therefore we can't render them until after the CSS transitions have happened.
 *    It also helps when we change the size of them (e.g. expand) them to simply get
 *    them to reinitialise with new dimensions
 */
const KronaPanel = ({config, viewOptions, canExpand, socket}) => {

    /* -----------    STATE MANAGEMENT    ------------------- */
    const [expanded, setExpanded] = useState(false);
    const [transitionInProgress, setTransitionInProgress] = useState(false);
    const transitionStarted = (duration=600) => { /* CSS transition is 0.5s */
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
    }

    /* ---------------   WHAT CHARTS DO WE RENDER?   -------------- */
    const renderCharts = () => {
        if (!expanded) return null;
        return (
            <div className="graphContainer">
                <div className="panelFlexRow">
                    <div krona_include="krona.html"></div>
                </div>
            </div>
        );
    }

    /* ----------------- R E N D E R ---------------- */
    const sampleColour = 'rgb(220,220,220)';
    return (
        <div
            className={`panelContainer ${expanded ? "expanded" : "collapsed"}`}
            style={{borderColor: sampleColour}}
        >
            <InfoRowKrona
                sampleColour={sampleColour}
                enableUserInteraction={canExpand}
                menuItems={menuItems}
                handleClick={toggleExpanded}
                isExpanded={expanded}
            />
            {transitionInProgress ? null : renderCharts()}
        </div>
    );
}

export default KronaPanel;

/*
 * Copyright (c) 2019 ARTIC Network http://artic.network
 * https://github.com/artic-network/rampart
 *
 * This file is part of RAMPART. RAMPART is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version. RAMPART is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
 * along with RAMPART. If not, see <http://www.gnu.org/licenses/>.
 *
 */

import React from 'react';
import {mouse, select} from "d3-selection";
import {calcScales} from "../../utils/commonFunctions";
import {heatColourScale} from "../../utils/colours";
import {getRelativeReferenceMapping} from "../../utils/config";
import Container, {Title, HoverInfoBox} from "./styles";

const EMPTY_CELL_COLOUR = "rgba(256, 256, 256, 0.15)"

/* given the DOM dimensions of the chart container, calculate the chart geometry (used by the SVG & D3) */
const calcChartGeom = (DOMRect) => ({
    width: DOMRect.width,
    height: DOMRect.height, // title line
    spaceLeft: DOMRect.width>600 ? 250 : // space for the reference names
        DOMRect.width>400 ? 150 :
            100,
    spaceRight: 0,
    spaceBottom: 90,
    spaceTop: 10,
    legendPadding: 20, /* horizontal */
    maxLegendWidth: 400
});

const calcCellDims = (chartGeom, numSamples, numReferences) => {
    const cellPadding = 1;
    const availableWidth = chartGeom.width - chartGeom.spaceLeft - chartGeom.spaceRight;
    const availableHeight = chartGeom.height - chartGeom.spaceBottom - chartGeom.spaceTop;
    const cellWidth = availableWidth/numSamples - cellPadding;
    const cellHeight = availableHeight/numReferences - cellPadding;
    return {
        height: cellHeight,
        width: cellWidth,
        padding: cellPadding
    }
}


const drawHeatMap = ({names, referencePanel, data, svg, scales, cellDims, chartGeom, relativeReferenceMapping, infoRef}) => {
    /* convert the refMatchPerSample data from raw counts to percentages & change to a d3-friendly struct.
    Input format:
      refMatchPerSample[sampleIdx][reference_idx] = INT
    Output data format:
      flat list, with each value itself a list:
        [sampleIdx, refPanelMatchIdx, fracIdentity]
    */
    console.log(`drawheatmap`);
    console.log(`with lengths ${names.length}`);
    console.log(`and ${referencePanel.length}`);
    const d3data = Array.from(new Array(names.length*referencePanel.length));

    let dataIdx = 0;

    let maxCount = 0;
    let total = 0;
    for (let sampleIdx=0; sampleIdx<names.length; sampleIdx++) {
        for (let refIdx=0; refIdx<referencePanel.length; refIdx++) {
            const count = parseInt(data[names[sampleIdx]].refMatches[referencePanel[refIdx].name]) || 0;
            if (count > maxCount) {
                maxCount = count;
            }
            total += count;
        }
    }

    // relativeReferenceMapping
    // if true then the heat is relative to the largest value, if false then it is the percentage
    // of reads by sample

    const showLegend = false;

    for (let sampleIdx=0; sampleIdx<names.length; sampleIdx++) {
        for (let refIdx=0; refIdx<referencePanel.length; refIdx++) {
            const count = parseInt(data[names[sampleIdx]].refMatches[referencePanel[refIdx].name]) || 0;
            const sampleTotal = parseInt(data[names[sampleIdx]].refMatches['total']) || 1;
            const percentOfSample = (100.0 * count) / sampleTotal;
            const percentOfTotal = (100.0 * count) / total;
            const heat = (100.0 * count) / (relativeReferenceMapping ? maxCount : sampleTotal);
            d3data[dataIdx] = {
                sampleIdx,
                refIdx,
                count,
                percentOfSample,
                percentOfTotal,
                heat
            };
            // console.log(names[sampleIdx] + " vs. " + referencePanel[refIdx].name + ": " + data[names[sampleIdx]].refMatches[referencePanel[refIdx].name] + " / " + data[names[sampleIdx]].refMatches['total'])
            dataIdx++;
        }
    }
    // /* NOTE scales.x(0) returns the far left pixel value of the cells, not the labels */

    /* remove the previous renderings... */

    svg.selectAll("*").remove();

    const referenceName = (d) => {
        const charPx = 8; /* guesstimate of character pixel width */
        const allowedChars = Math.floor(chartGeom.spaceLeft / charPx);
        if (d.name.length > allowedChars) {
            return `${d.name.slice(0,allowedChars-2)}...`;
        }
        return d.name;
    };

    /* render the reference names (on the far left) */
    svg.selectAll(".refLabel")
        .data(referencePanel) /* get the labels */
        .enter()
        .append("text")
        .attr("class", "refLabel axis")
        .text(referenceName)
        .attr('y', (refName, refIdx) => scales.y(refIdx+1) + 0.5*cellDims.height)
        .attr('x', chartGeom.spaceLeft - 8 /* - cellDims.height */)
        .attr("text-anchor", "end")
        .attr("font-size", "12px")
        .attr("alignment-baseline", "middle"); /* i.e. y value specifies top of text */

    // svg.selectAll(".refColour")
    //     .data(references) /* get the labels */
    //     .enter()
    //     .append("rect")
    //     .attr("class", "refColour")
    //     .attr('width', cellDims.height)
    //     .attr('height', cellDims.height)
    //     .attr("x", chartGeom.spaceLeft - 4 - cellDims.height)
    //     .attr('y', (refName, refIdx) => scales.y(refIdx+1))
    //     .attr("fill", (refName, refIdx) => referenceDiscreteColours[refIdx]);

    if (!showLegend) {
        const sampleName = (d) => {
            const charPx = 4; /* guesstimate of character pixel width */
            const allowedChars = Math.floor(chartGeom.spaceBottom / charPx);
            if (d.length > allowedChars) {
                return `${d.slice(0,allowedChars-2)}...`;
            }
            return d;
        };
        /* render the column labels (barcodes) on the bottom */
        svg.selectAll(".sampleName")
            .data(names)
            .enter()
            .append("g")
            .attr("class", "sampleNames axis")
            .append("text")
            .attr("class", "axis")
            .attr("transform", (name, idx) => `rotate(-45 ${scales.x(idx) + 0.5 * cellDims.width} ${chartGeom.height - chartGeom.spaceBottom + 10})`)
            .text(sampleName)
            .attr('x', (name, idx) => scales.x(idx) + 0.5 * cellDims.width)
            .attr('y', chartGeom.height - chartGeom.spaceBottom + 10)
            .attr("text-anchor", "end")
            .attr("font-size", "12px")
            .attr("alignment-baseline", "hanging");
    }

    function handleMouseMove(d, i) {
        const [mouseX, mouseY] = mouse(this); // [x, y] x starts from left, y starts from top
        const left  = mouseX > 0.5 * scales.x.range()[1] ? "" : `${mouseX + 16}px`;
        const right = mouseX > 0.5 * scales.x.range()[1] ? `${scales.x.range()[1] - mouseX}px` : "";
        // const mapString = referencePanel[d.refIdx].name !== "unmapped" ?
        //     `map to ${referencePanel[d.refIdx].name}` : `were not mapped to any reference`;
        const mapString = referencePanel[d.refIdx].name !== "unmapped" ?
            `Reference: ${referencePanel[d.refIdx].name}` : `Unmapped`;
        select(infoRef)
            .style("left", left)
            .style("right", right)
            .style("top", `${mouseY}px`)
            .style("visibility", "visible")
            .html(`
                Sample: ${names[d.sampleIdx]}
                <br/>
                ${mapString}
                <br/>
                ${d.count} reads
                <br/>
                ${d.percentOfSample.toFixed(2)}% of the sample
                <br/>
                ${d.percentOfTotal.toFixed(2)}% of the total
            `);
    }
    function handleMouseOut() {
        select(infoRef).style("visibility", "hidden");
    }

    /* render the coloured cells of the heatmap */
    svg.selectAll(".heatCell")
        .data(d3data)
        .enter()
        .append("rect")
        .attr("class", "heatCell")
        .attr('width', cellDims.width)
        .attr('height', cellDims.height)
        .attr("x", d => scales.x(d.sampleIdx) + cellDims.padding)
        .attr("y", d => scales.y(d.refIdx+1) + cellDims.padding)
        .attr("fill", d => d.count === 0 ? EMPTY_CELL_COLOUR : heatColourScale(d.heat))
        .on("mouseout", handleMouseOut)
        .on("mousemove", handleMouseMove);

    if (showLegend && relativeReferenceMapping) {
        /* render the legend (bottom) -- includes coloured cells & text */
        const legendDataValues = [0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
        let legendWidth = chartGeom.width - chartGeom.spaceRight - 2 * chartGeom.legendPadding;
        if (legendWidth > chartGeom.maxLegendWidth) legendWidth = chartGeom.maxLegendWidth;
        const legendRightOffset = (chartGeom.width - legendWidth) / 2;
        const legendBoxWidth = legendWidth / legendDataValues.length;
        const legendBoxHeight = 12;
        const legendRoof = chartGeom.height - chartGeom.spaceBottom + 30;
        const legendTextFn = (d, i) => {
            if (legendWidth === chartGeom.maxLegendWidth) return `${d}%`;
            if (i % 2) return `${d}%`;
            return "";
        };
        const legend = svg.selectAll(".legend")
            .data(legendDataValues)
            // .data(legendDataValues.slice(0, legendDataValues.length-1)) /* don't include the last one... */
            .enter().append("g")
            .attr("class", "legend");
        legend.append("rect")
            .attr('y', legendRoof)
            .attr("x", (d, i) => legendRightOffset + legendBoxWidth * (i + 1))
            .attr("width", legendBoxWidth)
            .attr("height", legendBoxHeight)
            .style("fill", (d) => {
                return d === 0 ? EMPTY_CELL_COLOUR : heatColourScale(d);
            });
        legend.append("text")
            .text(legendTextFn)
            .attr("class", "axis")
            .attr('x', (d, i) => legendRightOffset + legendBoxWidth * (i + 1))
            .attr('y', legendRoof + legendBoxHeight + 11)
            .attr("text-anchor", "middle")
            .attr("font-size", "12px")
            .attr("alignment-baseline", "hanging")
    }
};

class ReferenceHeatmap extends React.Component {
    constructor(props) {
        super(props);
        this.state = {chartGeom: {}, relativeReferenceMapping: false, hoverWidth: 0, svg: undefined};
    }
    redraw() {
        /* currently redo everything, but we could make this much much smarter */
        const sampleNames = Object.keys(this.props.data).filter((name) => name!=="all");
        const referencePanel = this.props.referencePanel.filter((info) => info.display);
        const chartGeom = this.state.chartGeom;
        const cellDims = calcCellDims(chartGeom, sampleNames.length, referencePanel.length);
        const scales = calcScales(
            chartGeom,
            sampleNames.length,     // number of columns
            referencePanel.length   // number of rows
        );
        drawHeatMap({
            names: sampleNames,
            referencePanel,
            data: this.props.data,
            svg: this.state.svg,
            scales,
            cellDims,
            chartGeom,
            relativeReferenceMapping: getRelativeReferenceMapping(this.props.config),
            infoRef: this.infoRef
        });
    }
    componentDidMount() {
        const svg = select(this.DOMref);
        const chartGeom = calcChartGeom(this.boundingDOMref.getBoundingClientRect());
        const hoverWidth = parseInt(chartGeom.width * 2/3, 10);
        this.setState({chartGeom, svg, hoverWidth})
    }
    componentDidUpdate() {
        this.redraw();
    }
    render() {
        return (
            <Container width={this.props.width} ref={(r) => {this.boundingDOMref = r}}>
                <Title>{this.props.title}</Title>
                <HoverInfoBox width={this.state.hoverWidth} ref={(r) => {this.infoRef = r}}/>
                <svg
                    ref={(r) => {this.DOMref = r}}
                    height={this.state.chartGeom.height || 0}
                    width={this.state.chartGeom.width || 0}
                />
                {this.props.renderProp ? this.props.renderProp : null}
            </Container>
        )
    }
}

export { ReferenceHeatmap };

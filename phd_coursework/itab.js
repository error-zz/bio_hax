

var row_ids;
var margin = { top: 80, right: 10, bottom: 50, left: 100 },
  cellSize=12;
  col_number=84;
  row_number=192;
  var row_ids;
  var selected_row_ids;
  width = cellSize*col_number, // - margin.left - margin.right,
  height = cellSize*row_number , // - margin.top - margin.bottom,
  //gridSize = Math.floor(width / 24),
  legendElementWidth = cellSize*2.5,
  colorBuckets = 4,
//  colors = ['#005824','#FFFFFF','#FFFF99','#91003F']; 

  colors = ['#336666','#FFFFFF','#FFEE66','#FF4444'];

  hcrow = [92, 53, 3, 16, 148, 13, 115, 105, 99, 97, 63, 48, 27, 25, 2, 19, 183, 180, 167, 142, 130, 95, 7, 6, 54, 35, 178, 175, 160, 158, 153, 145, 136, 133, 122, 118, 114, 112, 110, 108, 1, 82, 80, 79, 66, 45, 177, 172, 171, 166, 147, 141, 127, 120, 117, 102, 93, 46, 29, 18, 165, 159, 154, 146, 107, 184, 89, 83, 57, 149, 14, 132, 40, 31, 26, 134, 11, 77, 60, 56, 173, 185, 74, 62, 36, 34, 140, 121, 109, 85, 78, 174, 156, 150, 137, 104, 49, 138, 126, 10, 96, 76, 73, 38, 162, 51, 169, 168, 155, 15, 128, 5, 37, 20, 164, 94, 182, 181, 157, 131, 106, 72, 69, 55, 12, 64, 41, 143, 124, 59, 151, 144, 100, 67, 111, 75, 39, 186, 135, 43, 70, 68, 28, 188, 88, 4, 116, 187, 9, 24, 152, 33, 129, 113, 42, 139, 23, 103, 87, 61, 170, 90, 191, 161, 50, 123, 101, 32, 22, 17, 125, 71, 65, 91, 21, 192, 81, 52, 190, 179, 86, 47, 189, 84, 44, 163, 58, 119, 98, 8, 30, 176], 
  hccol = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84], 
 /* rowLabel = ['344', '299', '182', '260', '403', '257', '368', '357', '351', '349', '315', '294', '273', '271', '181', '264', '438', '435', '422', '397', '384', '347', '250', '249', '300', '281', '433', '430', '415', '413', '408', '400', '391', '387', '376', '372', '367', '364', '362', '360', '180', '334', '332', '331', '318', '291', '432', '427', '426', '421', '402', '396', '381', '374', '371', '354', '345', '292', '275', '263', '420', '414', '409', '401', '359', '439', '341', '335', '307', '404', '258', '386', '286', '277', '272', '388', '255', '329', '312', '306', '428', '440', '326', '314', '282', '280', '395', '375', '361', '337', '330', '429', '411', '405', '392', '356', '295', '393', '380', '254', '348', '328', '325', '284', '417', '297', '424', '423', '410', '259', '382', '248', '283', '266', '419', '346', '437', '436', '412', '385', '358', '324', '321', '305', '256', '316', '287', '398', '378', '311', '406', '399', '352', '319', '363', '327', '285', '441', '390', '289', '322', '320', '274', '443', '340', '210', '369', '442', '252', '270', '407', '279', '383', '366', '288', '394', '269', '355', '339', '313', '425', '342', '446', '416', '296', '377', '353', '278', '268', '262', '379', '323', '317', '343', '267', '447', '333', '298', '445', '434', '338', '293', '444', '336', '290', '418', '310', '373', '350', '251', '276', '431'], // change to gene name or probe id */
  rowLabel = ['222-rdj', '222-r6k', '222-kcp', '222-n7w', '222-w6s', '222-njd', '222-trb', '222-scr', '222-tm2', '222-rdm', '222-s5t', '222-qdg', '222-nxs', '222-pd2', '222-k4w', '222-nn8', '223-25n', '223-2dq', '222-xjs', '222-w64', '222-vtc', '222-qr7', '222-kqc', '222-kpb', '222-q85', '222-ppc', '222-zt7', '222-vxk', '222-xg6', '222-xg7', '222-w53', '222-wt9', '222-wcg', '222-vcp', '222-vw5', '222-t3f', '222-v6v', '222-tcs', '222-sqg', '222-txd', '222-kcn', '222-rf3', '222-sbt', '222-r5p', '222-s8r', '222-pzb', '222-zg9', '222-xxs', '222-xgz', '222-xrq', '222-wqj', '222-w5w', '222-vxm', '222-vpf', '222-s5f', '222-shc', '222-tbv', '222-phx', '222-q6q', '222-nk9', '222-xkm', '222-x9j', '222-w8k', '222-ww5', '222-v2w', '222-2xf', '222-t76', '222-sqn', '222-rvh', '222-wd6', '222-mtv', '222-w5j', '222-pkd', '222-q5s', '222-nt9', '222-w6w', '222-n54', '222-rbm', '222-qmb', '222-qg3', '222-xxn', '233-2kk', '222-scf', '222-qs4', '222-q86', '222-pfn', '222-vxt', '222-vqm', '222-tpm', '222-sch', '222-spj', '222-z2j', '222-xc2', '222-x2p', '222-w9n', '222-smn', '222-ppm', '222-v85', '222-v75', '222-mzm', '222-thk', '222-qz4', '222-scq', '222-pz9', '222-xgr', '222-qsw', '222-z2p', '222-z34', '222-x97', '222-nk8', '222-w4f', '222-266', '222-pvz', '222-nxg', '222-xnv', '222-tb2', '222-z9v', '223-2cn', '222-x7v', '222-w4d', '222-v2k', '222-r5n', '222-q98', '222-qr5', '222-mgg', '222-s5g', '222-pgx', '222-wmc', '222-t3x', '222-ppj', '222-x4v', '222-wv7', '222-smp', '222-s9z', '222-trc', '222-smt', '222-qs6', '223-35n', '222-w9q', '222-qvg', '222-s9h', '222-r44', '222-pz8', '223-45g', '222-t8p', '222-kbj', '222-tnk', '223-2r9', '222-mts', '222-pq7', '222-x43', '222-qfh', '222-w58', '222-t3m', '222-qs8', '222-w6m', '222-p93', '222-tq2', '222-s9x', '222-qst', '222-xgx', '222-scd', '223-3zs', '222-z4d', '222-qsv', '222-vwn', '222-trk', '222-pkc', '222-pmc', '222-nsr', '222-tsx', '222-sbh', '222-qz5', '222-sxc', '222-pc3', '223-3hw', '222-s5k', '222-ppf', '223-5fx', '222-zxc', '222-t4z', '222-r56', '223-42k', '222-sqm', '222-qz8', '222-xpx', '222-rx9', '222-vpc', '222-tk9', '222-mr7', '222-q7v', '222-z64'];
  colLabel = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84']; // change to contrast name

d3.tsv("devel_data_heatmap.tsv",
function(d) {
  return {
    row:   +d.row_idx,
    col:   +d.col_idx,
    value: +d.log2ratio
  };
},
function(error, data) {
  var colorScale = d3.scale.quantile()
      .domain([ -10 , 0, 2, 10])
      .range(colors);
  
  var svg = d3.select("#chart").append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      ;
  var rowSortOrder=false;
  var colSortOrder=false;
  var rowLabels = svg.append("g")
      .selectAll(".rowLabelg")
      .data(rowLabel)
      .enter()
      .append("text")
      .text(function (d) { return d; })
      .attr("x", 0)
      .attr("y", function (d, i) { return hcrow.indexOf(i+1) * cellSize; })
      .style("text-anchor", "end")
      .attr("transform", "translate(-6," + cellSize / 1.5 + ")")
      .attr("class", function (d,i) { return "rowLabel mono r"+i;} ) 
      .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
      .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
      .on("click", function(d,i) {rowSortOrder=!rowSortOrder; sortbylabel("r",i,rowSortOrder);d3.select("#order").property("selectedIndex", 4).node().focus();;})
      ;

  var colLabels = svg.append("g")
      .selectAll(".colLabelg")
      .data(colLabel)
      .enter()
      .append("text")
      .text(function (d) { return d; })
      .attr("x", 0)
      .attr("y", function (d, i) { return hccol.indexOf(i+1) * cellSize; })
      .style("text-anchor", "left")
      .attr("transform", "translate("+cellSize/2 + ",-6) rotate (-90)")
      .attr("class",  function (d,i) { return "colLabel mono c"+i;} )
      .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
      .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
      .on("click", function(d,i) {colSortOrder=!colSortOrder;  sortbylabel("c",i,colSortOrder);d3.select("#order").property("selectedIndex", 4).node().focus();;})
      ;

  var heatMap = svg.append("g").attr("class","g3")
        .selectAll(".cellg")
        .data(data,function(d){return d.row+":"+d.col;})
        .enter()
        .append("rect")
        .attr("x", function(d) { return hccol.indexOf(d.col) * cellSize; })
        .attr("y", function(d) { return hcrow.indexOf(d.row) * cellSize; })
        .attr("class", function(d){return "cell cell-border cr"+(d.row-1)+" cc"+(d.col-1);})
        .attr("width", cellSize)
        .attr("height", cellSize)
        .style("fill", function(d) { return colorScale(d.value); })
        /* .on("click", function(d) {
               var rowtext=d3.select(".r"+(d.row-1));
               if(rowtext.classed("text-selected")==false){
                   rowtext.classed("text-selected",true);
               }else{
                   rowtext.classed("text-selected",false);
               }
        })*/
        .on("mouseover", function(d){
               //highlight text
               d3.select(this).classed("cell-hover",true);
               d3.selectAll(".rowLabel").classed("text-highlight",function(r,ri){ return ri==(d.row-1);});
               d3.selectAll(".colLabel").classed("text-highlight",function(c,ci){ return ci==(d.col-1);});
        
               //Update the tooltip position and value
               d3.select("#tooltip")
                 .style("left", (d3.event.pageX+10) + "px")
                 .style("top", (d3.event.pageY-10) + "px")
                 .select("#value")
                 .text("lables:"+rowLabel[d.row-1]+","+colLabel[d.col-1]+"\ndata:"+d.value+"\nrow-col-idx:"+d.col+","+d.row+"\ncell-xy "+this.x.baseVal.value+", "+this.y.baseVal.value);  
               //Show the tooltip
               d3.select("#tooltip").classed("hidden", false);
        })
        .on("mouseout", function(){
               d3.select(this).classed("cell-hover",false);
               d3.selectAll(".rowLabel").classed("text-highlight",false);
               d3.selectAll(".colLabel").classed("text-highlight",false);
               d3.select("#tooltip").classed("hidden", true);
                            
                //var qqq_selected = d3.select(this).classed("text-hover",true) || null;
                // document.write(qqq_selected);
             // var to = date.top(1)[0].visit_date.format("YYYY-MM-DD") || null;
             // headers = [title, "Subset Count <br><small> " + from + " to " + to + "</small", "Subset Percentage <br> <small>" + from + " to " + to + "</small>", "Total Counts", "Total Pecentage"]
  
  
        })
        ;

  var legend = svg.selectAll(".legend")
      .data([-10,0,2,10])
      .enter().append("g")
      .attr("class", "legend");
 
  legend.append("rect")
    .attr("x", function(d, i) { return legendElementWidth * i; })
    .attr("y", height+(cellSize*2))
    .attr("width", legendElementWidth)
    .attr("height", cellSize)
    .style("fill", function(d, i) { return colors[i]; });
 
  legend.append("text")
    .attr("class", "mono")
    .text(function(d) { return d; })
    .attr("width", legendElementWidth)
    .attr("x", function(d, i) { return legendElementWidth * i; })
    .attr("y", height + (cellSize*4));

// Change ordering of cells

  function sortbylabel(rORc,i,sortOrder){
       var t = svg.transition().duration(1000);
       var log2r=[];
       var sorted; // sorted is zero-based index
       d3.selectAll(".c"+rORc+i) 
         .filter(function(ce){
            log2r.push(ce.value);
          })
       ;
       if(rORc=="r"){ // sort log2ratio of a gene
         sorted=d3.range(col_number).sort(function(a,b){ if(sortOrder){ return log2r[b]-log2r[a];}else{ return log2r[a]-log2r[b];}});
         t.selectAll(".cell")
           .attr("x", function(d) { return sorted.indexOf(d.col-1) * cellSize; })
           ;
         t.selectAll(".colLabel")
          .attr("y", function (d, i) { return sorted.indexOf(i) * cellSize; })
         ;
       }else{ // sort log2ratio of a contrast
         sorted=d3.range(row_number).sort(function(a,b){if(sortOrder){ return log2r[b]-log2r[a];}else{ return log2r[a]-log2r[b];}});
         t.selectAll(".cell")
           .attr("y", function(d) { return sorted.indexOf(d.row-1) * cellSize; })
           ;
         t.selectAll(".rowLabel")
          .attr("y", function (d, i) { return sorted.indexOf(i) * cellSize; })
         ;
       }
  }

  d3.select("#order").on("change",function(){
    order(this.value);
  });
  
  function order(value){
   if(value=="hclust"){
    var t = svg.transition().duration(1500);
    t.selectAll(".cell")
      .attr("x", function(d) { return hccol.indexOf(d.col) * cellSize; })
      .attr("y", function(d) { return hcrow.indexOf(d.row) * cellSize; })
      ;

    t.selectAll(".rowLabel")
      .attr("y", function (d, i) { return hcrow.indexOf(i+1) * cellSize; })
      ;

    t.selectAll(".colLabel")
      .attr("y", function (d, i) { return hccol.indexOf(i+1) * cellSize; })
      ;

   }else if (value=="probecontrast"){
    var t = svg.transition().duration(1500);
    t.selectAll(".cell")
      .attr("x", function(d) { return (d.col - 1) * cellSize; })
      .attr("y", function(d) { return (d.row - 1) * cellSize; })
      ;

    t.selectAll(".rowLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;

    t.selectAll(".colLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;

   }else if (value=="probe"){
    var t = svg.transition().duration(1500);
    t.selectAll(".cell")
      .attr("y", function(d) { return (d.row - 1) * cellSize; })
      ;

    t.selectAll(".rowLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;
   }else if (value=="contrast"){
    var t = svg.transition().duration(1500);
    t.selectAll(".cell")
      .attr("x", function(d) { return (d.col - 1) * cellSize; })
      ;
    t.selectAll(".colLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;
   }
  }
  
  var qqq; 
  var sa=d3.select(".g3")
      .on("mousedown", function() {
          if( !d3.event.altKey) {
             d3.selectAll(".cell-selected").classed("cell-selected",false);
             d3.selectAll(".rowLabel").classed("text-selected",false);
             d3.selectAll(".colLabel").classed("text-selected",false);
          }
         var p = d3.mouse(this);
         sa.append("rect")
         .attr({
             rx      : 0,
             ry      : 0,
             class   : "selection",
             x       : p[0],
             y       : p[1],
             width   : 1,
             height  : 1
         })
      })
      .on("mousemove", function() {
         var s = sa.select("rect.selection");
      
         if(!s.empty()) {
             var p = d3.mouse(this),
                 d = {
                     x       : parseInt(s.attr("x"), 10),
                     y       : parseInt(s.attr("y"), 10),
                     width   : parseInt(s.attr("width"), 10),
                     height  : parseInt(s.attr("height"), 10)
                 },
                 move = {
                     x : p[0] - d.x,
                     y : p[1] - d.y
                 }
             ;
             
             // JM - 5/6 - EDIT
             // document.write('<p>x: ' + d.x + '<br>y: ' + d.y + '<br>width: ' + d.width + '<br>height: ' + d.height + '</p>');
             function dump(){
                    document.getElementById('d_x').innerHTML = d.x;
                    document.getElementById('d_y').innerHTML = d.y;
                    document.getElementById('d_width').innerHTML = d.width;
                    document.getElementById('d_height').innerHTML = d.height;
                    // document.getElementById('d_rows').innerHTML =  s;
                    // document.getElementById('row_ids').innerHTML = d3.selectAll(".rowLabel.text-selected");
                    // document.getElementById('col_ids').innerHTML = d3.selectAll(".colLabel.text-selected");
                    // console.log("TEST");
                    // console.log(d3.selectAll(".rowLabel.text-selected")[0].map(function(e){
                        // console.log(e.textContent);
                        // return e.textContent;
                    // }));
                    row_ids = d3.selectAll(".rowLabel.text-selected")[0].map(function(e){
                         //console.log(e.textContent);
                         return e.textContent;
                    });
                    document.getElementById('row_ids').innerHTML = row_ids.join("<br>");
                    document.getElementById('row_ids2').innerHTML = row_ids.join("<br>");
                    console.log(row_ids);

                    var col_ids = d3.selectAll(".colLabel.text-selected")[0].map(function(e){
                         //console.log(e.textContent);
                         return e.textContent;
                    });
                    document.getElementById('col_ids').innerHTML = col_ids.join("<br>");
                    document.getElementById('col_ids2').innerHTML = col_ids.join("<br>");

                    // document.getElementById('col_ids').innerHTML = d3.selectAll(".colLabel").classed("text-selected",true);
                    // document.getElementById('col_ids').innerHTML = d3.selectAll(".colLabel").classed("text-selected",false);
                    // s.classed("cell-selected", true);
            }
            dump();
            //  JM - 5/6 - EDIT
      
             if(move.x < 1 || (move.x*2<d.width)) {
                 d.x = p[0];
                 d.width -= move.x;
             } else {
                 d.width = move.x;       
             }
      
             if(move.y < 1 || (move.y*2<d.height)) {
                 d.y = p[1];
                 d.height -= move.y;
             } else {
                 d.height = move.y;       
             }
             s.attr(d);
      
                 // deselect all temporary selected state objects
             d3.selectAll('.cell-selection.cell-selected').classed("cell-selected", false);
             d3.selectAll(".text-selection.text-selected").classed("text-selected",false);

             d3.selectAll('.cell').filter(function(cell_d, i) {
                 if(
                     !d3.select(this).classed("cell-selected") && 
                         // inner circle inside selection frame
                     (this.x.baseVal.value)+cellSize >= d.x && (this.x.baseVal.value)<=d.x+d.width && 
                     (this.y.baseVal.value)+cellSize >= d.y && (this.y.baseVal.value)<=d.y+d.height
                 ) {
      
                     d3.select(this)
                     .classed("cell-selection", true)
                     .classed("cell-selected", true);

                     d3.select(".r"+(cell_d.row-1))
                     .classed("text-selection",true)
                     .classed("text-selected",true);

                     d3.select(".c"+(cell_d.col-1))
                     .classed("text-selection",true)
                     .classed("text-selected",true);
                 }
             });
         }
      })
      .on("mouseup", function() {
            // remove selection frame
         sa.selectAll("rect.selection").remove();
      
             // remove temporary selection marker class
         d3.selectAll('.cell-selection').classed("cell-selection", false);
         d3.selectAll(".text-selection").classed("text-selection",false);
      })
      .on("mouseout", function() {
         if(d3.event.relatedTarget.tagName=='html') {
                 // remove selection frame
             sa.selectAll("rect.selection").remove();
                 // remove temporary selection marker class
             d3.selectAll('.cell-selection').classed("cell-selection", false);
             d3.selectAll(".rowLabel").classed("text-selected",false);
             d3.selectAll(".colLabel").classed("text-selected",false);
         }
      })
      ;
});


/* PREVIOUSLY CONTENT OF DEMOGRAPHICS.JS */

function createDemographicFilters(codebook) {

  // Get the demographic values and choice lists
  function filterDemographicValues(element) {
    return element.table == "PatientRegistrationAndDemographics" && element.type == "choice";
  }

  var demographic_codes = codebook.filter(filterDemographicValues);
  return demographic_codes;

}

function getRaceMap(codebook, field) {
  // Fetch race code
  var race_map = codebook.filter(function(d) { return d.field == field});

  keys = {};

  race_map[0].choices.split(";").map(function(d) {
    d.split("=").reduce( function(p,d) { keys[p] = d });
  });
  
  // Explode and create keys
  return keys;
}

queue()
.defer(d3.csv, "static/cctg/PatientRegistrationAndDemographics.csv")
.defer(d3.csv, "static/cctg/visit.csv")
.defer(d3.csv, "static/cctg/codebook.csv")
.await(function(error, demographics, visits, codebook) { 

  // Lambda functions
  function isEnrollmentAnd595(element) {
    return element.cycle == 0 && element.study_name == "595";
  }

  function adaptVisitDates(element, i) {
    element.visit_date = moment(element.visit_date);
    return element;
  }

  function addIndex(element, i) {
    element.index = i + 1;
    return element;
  }
  
  var demographic_codes = createDemographicFilters(codebook);
  var race_map = getRaceMap(demographic_codes, "race");
  var gender_map =  getRaceMap(demographic_codes, "gender");
  var ethnicity_map = getRaceMap(demographic_codes, "ethn");

  // Main
  var enrollments = visits.filter(isEnrollmentAnd595);
  var enrollments = enrollments.map(adaptVisitDates);


  // JM update

  // console.log(row_ids2);

  var sorted_enrollments = enrollments.sort(function (a, b) {
      if (a.visit_date.isAfter(b.visit_date))
        return 1;
      if (a.visit_date.isBefore(b.visit_date))
        return -1;
      return 0;
  });

  var sorted_enrollments = enrollments.map(addIndex);

  // We need to join demographic and enrollment data
  var joined_table = demographics.reduce(function(p,d) {
    p[d.pid] = d;
    return p;
  }, {});

  sorted_enrollments = sorted_enrollments.map(function(d) { 
      if (d.our in joined_table) {
        d.ethn = joined_table[d.our].ethn || "";
        d.race = joined_table[d.our].race || "";
        d.gender = joined_table[d.our].gender || "";
      }
      return d;
  })

  // Create enrollments crossfilter
  var cf_enrollment = crossfilter(sorted_enrollments),
      all = cf_enrollment.groupAll(),
      date = cf_enrollment.dimension(function(d) { return d.visit_date.toDate(); }),
      dates = date.group(d3.time.day);
      ethnicity = cf_enrollment.dimension(function(d) { return d.ethn; }),
      ethnicities = ethnicity.group(),
      race = cf_enrollment.dimension(function(d) { return d.race; }),
      races = race.group(),
      gender = cf_enrollment.dimension(function(d) { return d.gender; }),
      genders = gender.group(),
      total_race = cf_enrollment.dimension(function(d) { return d.race; }),
      total_races = total_race.group(),
      total_gender = cf_enrollment.dimension(function(d) { return d.gender; }),
      total_genders = total_gender.group(),
      total_ethn = cf_enrollment.dimension(function(d) { return d.ethn; }),
      total_ethns = total_ethn.group();

  const ethn_totals = total_ethns.top(Infinity);
  const ethn_count_total = total_ethn.groupAll().reduceCount().value();
  total_ethns.dispose();

  const gender_totals = total_genders.top(Infinity);
  const gender_count_total = total_gender.groupAll().reduceCount().value();
  total_genders.dispose();

  const race_totals = total_races.top(Infinity);
  const race_count_total = total_race.groupAll().reduceCount().value();
  total_races.dispose();

  var table = d3.selectAll(".table")
      .data([function(d) {return demographicsList(d, date, "Race", race, race_map, race_totals, race_count_total)}]);

  $("#race-btn").bind("click",function(){
    $(".btn").removeClass("active");
    $(this).addClass("active");
    var table = d3.selectAll(".table")
        .data([function(d) {return demographicsList(d, date, "Race", race, race_map, race_totals, race_count_total)}]);
    renderAll()  
  });


  $("#ethn-btn").bind("click",function(){
    $(".btn").removeClass("active");
    $(this).addClass("active");
    var table = d3.selectAll(".table")
        .data([function(d) {return demographicsList(d, date, "Ethnicity", ethnicity, ethnicity_map, ethn_totals, ethn_count_total)}]);
    renderAll()  
  });

  $("#gender-btn").bind("click",function(){
    $(".btn").removeClass("active");
    $(this).addClass("active");
    var table = d3.selectAll(".table")
        .data([function(d) {return demographicsList(d, date, "Gender", gender, gender_map, gender_totals, gender_count_total)}]);
    renderAll()  
  });



  var earliest = date.bottom(1)[0].visit_date.toDate();
  var latest = date.top(1)[0].visit_date.toDate();

  lineChart(sorted_enrollments, "#accrual-chart");

  var charts = [
    barChart()
        .dimension(date)
        .group(dates)
        .round(d3.time.day.round)
      .x(d3.time.scale()
        .domain([earliest,latest])
        .rangeRound([0, 10 * 90]))
        .filter([new Date(2013, 10, 1), new Date(2014, 4, 1)])
        
    // lineChart()
    //    .dimension(date)
    //    .group(dates)
    //    .round(d3.time.day.round)
    //  .x(d3.time.scale()
    //    .domain([earliest,latest])
    //    .rangeRound([0, 10 * 90]))
    //    .filter([new Date(2013, 10, 1), new Date(2014, 4, 1)])
        
  ];

  function renderAll() {
    chart.each(render);
    table.each(render);
  }

  // Renders the specified chart or list.
  function render(method) {
    d3.select(this).call(method);
  }

  var chart = d3.selectAll(".chart")
      .data(charts)
      .each(function(chart) { chart.on("brush", renderAll).on("brushend", renderAll); });

  renderAll();

});

function lineChart(data, element) {

  var margin = {top: 20, right: 20, bottom: 30, left: 50},
      width = 960 - margin.left - margin.right,
      height = 500 - margin.top - margin.bottom;

  var x = d3.time.scale()
      .range([0, width]);


  var brush = d3.svg.brush()
    .x(x)
    .y(y)
    .extent([new Date(2014, 1, 2), new Date(2014, 4, 3)])
    .on("brushend", brushended);

  var y = d3.scale.linear()
      .range([height, 0]);

  var xAxis = d3.svg.axis()
      .scale(x)
      .tickFormat(d3.time.format("%Y - %m"))
      .orient("bottom");

  var yAxis = d3.svg.axis()
      .scale(y)
      .orient("left");

  var line = d3.svg.line()
      .x(function(d) { return x(d.visit_date.toDate()); })
      .y(function(d) { return y(d.index); });

  var svg = d3.select(element).append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  x.domain(d3.extent(data, function(d) { return d.visit_date.toDate(); }));
  y.domain(d3.extent(data, function(d) { return d.index; }));

  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis)
    .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("Cumulative Accrual")

  //var gBrush = svg.append("g")
  //    .attr("class", "brush")
  //    .call(brush)
  //    .call(brush.event);

  //gBrush.selectAll("rect")
  //    .attr("height", height);

  svg.append("path")
      .datum(data)
      .attr("class", "line")
      .attr("d", line);

}

function brushed() {
  var extent = brush.extent();
  point.each(function(d) { d.selected = false; });
  search(quadtree, extent[0][0], extent[0][1], extent[1][0], extent[1][1]);
  point.classed("selected", function(d) { return d.selected; });
}

function brushended() {
  if (!d3.event.sourceEvent) return; // only transition after input
  var extent0 = brush.extent(),
      extent1 = extent0.map(d3.time.day.round);

  // if empty when rounded, use floor & ceil instead
  if (extent1[0] >= extent1[1]) {
    extent1[0] = d3.time.day.floor(extent0[0]);
    extent1[1] = d3.time.day.ceil(extent0[1]);
  }

  d3.select(this).transition()
      .call(brush.extent(extent1))
      .call(brush.event);
}

function barChart() {
  if (!barChart.id) barChart.id = 0;

  var margin = {top: 20, right: 10, bottom: 30, left: 40},
      x,
      y = d3.scale.linear().range([100, 0]),
      id = barChart.id++,
      axis = d3.svg.axis().orient("bottom").tickFormat(d3.time.format("%Y - %m")),
      brush = d3.svg.brush(),
      brushDirty,
      dimension,
      group,
      round;

  function chart(div) {
    var width = x.range()[1],
        height = y.range()[0];

    y.domain([0, group.top(1)[0].value]);

    div.each(function() {
      var div = d3.select(this),
          g = div.select("g");

      // Create the skeletal chart.
      if (g.empty()) {
        div.select(".title").append("a")
            .attr("href", "javascript:reset(" + id + ")")
            .attr("class", "reset")
            .text("reset")
            .style("display", "none");

        g = div.append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
          .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        g.append("clipPath")
            .attr("id", "clip-" + id)
          .append("rect")
            .attr("width", width)
            .attr("height", height);

        g.selectAll(".bar")
            .data(["background", "foreground"])
          .enter().append("path")
            .attr("class", function(d) { return d + " bar"; })
            .datum(group.all());

        g.selectAll(".foreground.bar")
            .attr("clip-path", "url(#clip-" + id + ")");

        g.append("g")
            .attr("class", "axis")
            .attr("transform", "translate(0," + height + ")")
            .call(axis);

        // Initialize the brush component with pretty resize handles.
        var gBrush = g.append("g").attr("class", "brush").call(brush);
        gBrush.selectAll("rect").attr("height", height);
        gBrush.selectAll(".resize").append("path").attr("d", resizePath);
      }

      // Only redraw the brush if set externally.
      if (brushDirty) {
        brushDirty = false;
        g.selectAll(".brush").call(brush);
        div.select(".title a").style("display", brush.empty() ? "none" : null);
        if (brush.empty()) {
          g.selectAll("#clip-" + id + " rect")
              .attr("x", 0)
              .attr("width", width);
        } else {
          var extent = brush.extent();
          g.selectAll("#clip-" + id + " rect")
              .attr("x", x(extent[0]))
              .attr("width", x(extent[1]) - x(extent[0]));
        }
      }

      g.selectAll(".bar").attr("d", barPath);
    });

    function barPath(groups) {
      var path = [],
          i = -1,
          n = groups.length,
          d;
      while (++i < n) {
        d = groups[i];
        path.push("M", x(d.key), ",", height, "V", y(d.value), "h9V", height);
      }
      return path.join("");
    }

    function resizePath(d) {
      var e = +(d == "e"),
          x = e ? 1 : -1,
          y = height / 3;
      return "M" + (.5 * x) + "," + y
          + "A6,6 0 0 " + e + " " + (6.5 * x) + "," + (y + 6)
          + "V" + (2 * y - 6)
          + "A6,6 0 0 " + e + " " + (.5 * x) + "," + (2 * y)
          + "Z"
          + "M" + (2.5 * x) + "," + (y + 8)
          + "V" + (2 * y - 8)
          + "M" + (4.5 * x) + "," + (y + 8)
          + "V" + (2 * y - 8);
    }
  }

  brush.on("brushstart.chart", function() {
    var div = d3.select(this.parentNode.parentNode.parentNode);
    div.select(".title a").style("display", null);
  });

  brush.on("brush.chart", function() {
    var g = d3.select(this.parentNode),
        extent = brush.extent();
    if (round) g.select(".brush")
        .call(brush.extent(extent = extent.map(round)))
      .selectAll(".resize")
        .style("display", null);
    g.select("#clip-" + id + " rect")
        .attr("x", x(extent[0]))
        .attr("width", x(extent[1]) - x(extent[0]));
    dimension.filterRange(extent);
  });

  brush.on("brushend.chart", function() {
    if (brush.empty()) {
      var div = d3.select(this.parentNode.parentNode.parentNode);
      div.select(".title a").style("display", "none");
      div.select("#clip-" + id + " rect").attr("x", null).attr("width", "100%");
      dimension.filterAll();
    }
  });

  chart.margin = function(_) {
    if (!arguments.length) return margin;
    margin = _;
    return chart;
  };

  chart.x = function(_) {
    if (!arguments.length) return x;
    x = _;
    axis.scale(x);
    brush.x(x);
    return chart;
  };

  chart.y = function(_) {
    if (!arguments.length) return y;
    y = _;
    return chart;
  };

  chart.dimension = function(_) {
    if (!arguments.length) return dimension;
    dimension = _;
    return chart;
  };

  chart.filter = function(_) {
    if (_) {
      brush.extent(_);
      dimension.filterRange(_);
    } else {
      brush.clear();
      dimension.filterAll();
    }
    brushDirty = true;
    return chart;
  };

  chart.group = function(_) {
    if (!arguments.length) return group;
    group = _;
    return chart;
  };

  chart.round = function(_) {
    if (!arguments.length) return round;
    round = _;
    return chart;
  };

  return d3.rebind(chart, brush, "on");
}

function demographicsList(table, date, title, dimension, codes, all_values, all_totals) {

  var dimension_group = dimension.group().reduceCount();
  //var demographics_values = dimension_group.top(Infinity);
  var demographics_values = dimension_group.all();
  var total = dimension.groupAll().reduceCount().value();

  // Make percentages
  demographics_values.forEach(function(d, i) {
    d.index = i + 1;
    d.percent = d.value/total;
  });

  // Join demographics values with total values
  demographics_values.forEach(function(d, i) {
    d.total_count = all_values.filter(function(f){return d.key == f.key})[0].value || 0;
    d.total_percent = d.total_count/all_totals;
  });

  var demographics_values = demographics_values.sort(function (a, b) {
      if (parseInt(a.key) > parseInt(b.key))
        return 1;
      if (parseInt(a.key) < parseInt(b.key))
        return -1;
      return 0;
  });


  // Labels
  var from = date.bottom(1)[0].visit_date.format("YYYY-MM-DD") || null;
  var to = date.top(1)[0].visit_date.format("YYYY-MM-DD") || null;
  headers = [title, "Subset Count <br><small> " + from + " to " + to + "</small", "Subset Percentage <br> <small>" + from + " to " + to + "</small>", "Total Counts", "Total Pecentage"]

  table.each(function() {

    var thead = d3.select(".table").selectAll("thead");

    var ths = thead.selectAll("th");
    ths.remove();

    var ths = thead.selectAll("th");
    var th = ths
        .data(headers)
        .enter().append("th");

    var ith = th.append("th")
          .html(function(d) { return d; });

    var tbody = d3.select(this).select("tbody");

    // Do not remove headers
    var trs = tbody.selectAll("tr");
    trs.remove();

    var trs = tbody.selectAll("tr");

    var tr = trs 
        .attr("class", "top-ten-item")
        .data(demographics_values)
        .enter().append("tr");

    var td = tr.append("td")
          .text(function(d) { 
              var all_keys = d.key.split(";");
              str = ""
              all_keys.forEach(function(p, i) {
                  str += codes[p] || "No Answer"
                  if(i < all_keys.length - 1) {
                    str += " and ";
                  }
              })
              return str;
            })
            .attr("class", "col-md-6");

    var td = tr.append("td")
          .text(function(d) { return d.value; })
          .attr("class", "col-md-2");

    var td = tr.append("td")
          .text(function(d) { return d3.format("%")(d.percent); })
          .attr("class", "col-md-2");

    var td = tr.append("td")
          .text(function(d) { return d.total_count; })
          .attr("class", "col-md-2");

    var td = tr.append("td")
          .text(function(d) { return d3.format("%")(d.total_percent); })
          .attr("class", "col-md-2");

    var sum = 0;
    demographics_values.map(function(d) { sum = sum + d.value});

    var per_sum = 0;
    demographics_values.map(function(d) { per_sum = per_sum + d.percent});

    var tfoot = d3.select(this).select("tfoot");

    var trs = tfoot.selectAll("tr");
    trs.remove();

    //var trs = tfoot.selectAll("tr");

    var tr = tfoot.append("tr");

    var td = tr.append("td")
          .html("<b>Sum</b>")
          .attr("class", "col-md-6");

    var td = tr.append("td")
          .text(sum)
          .attr("class", "col-md-2");

    var td = tr.append("td")
          .text(d3.format("%")(per_sum))
          .attr("class", "col-md-2");

    var td = tr.append("td")
          .text(all_totals)
          .attr("class", "col-md-2");

    var td = tr.append("td")
          .text(d3.format("%")(1))
          .attr("class", "col-md-2");
  });

}

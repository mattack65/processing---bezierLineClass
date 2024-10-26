// Author: matt@void.ch
// This class represents a C1 continuous, smooth Bezier curve through all the points that you give it
// The data contains the anchor-points (given by you) and the control-points (calculated by the class)

// Author: matt@void.ch
// This class represents a C1 continuous, smooth Bezier curve through all the points that you give it
// The data contains the anchor-points (given by you) and the control-points (calculated by the class)

// todos: 
// 1) make the swell factor changeable while the bezier line already exists
// 2) test and make move() more efficient
// 3) prevent poinst from being added after the curve has been closed or alternatively implement that case correctly
// 4) allow sthe curve to be opened again after it had been closed
// 5) allow the retrieval of other point arrays than just equidistant ones. Allow the handing over of a function that maps 0 ... 1 to 0 .. 1 and returns points in those mapped intervals
// 6) merging teo bezier curves?
// 7) find the intersectionpoints between two curves or even within one curve? needed?


import java.util.Locale;

// todo: we could make this now hard-coded factor a parameter and therefore change the "curviness" of the lines.
// it makes quite a massive difference to how the lines look

class BezierLine {
  ArrayList<PVector> anchorPoints;   // these are the anchor points of the curve
  ArrayList<PVector> ctrlPointsL;    // these are the control points to the left of each anchor point
  ArrayList<PVector> ctrlPointsR;    // these are the control points to the right of each anchor point
  
  boolean m_closedCurve;       // initially we assume that the curve is not closed
  float m_swell;               // swell describes how lose or how tight the curce winds around the anchor points. 
                               // the bigger the number the more oppulent the curves. 
                               // More than about 0.8 is porbably too much. 0.3 is default.  
  
  BezierLine() {
    this(0.3f);                // use the default swell
  }
    
  BezierLine(float swell) {
    anchorPoints = new ArrayList<PVector>();
    ctrlPointsL = new ArrayList<PVector>();
    ctrlPointsR = new ArrayList<PVector>();
    m_swell = swell;
    m_closedCurve = false;
  }

  void addPoint(PVector point) {
    // add a single anchor point 
    anchorPoints.add(point);
    updateControlPoints();
  }
  
  void movePoint(int index, PVector moveVector) {
    // move a single anchor point by a vector 
    anchorPoints.get(index).x += moveVector.x;
    anchorPoints.get(index).y += moveVector.y;
    // recalc everything (it's complete overkill: todo: only re-calc the controlpoints necessary
    ctrlPointsL.clear();
    ctrlPointsR.clear();
    updateControlPoints(); 
  }

  void addPoints(PVector[] newPoints) {
    // add an array of points 
    for (PVector p : newPoints) {
      anchorPoints.add(p);
    }
    updateControlPoints();  
  }
  
  void addPoints(ArrayList<PVector> newPointsList) {
    // add a list of points 
    for (PVector p : newPointsList) {
      anchorPoints.add(p);
    }
    updateControlPoints();  
  }
  
  void closeCurve() {
    if (m_closedCurve || anchorPoints.size() < 3) return; // Ensure there are enough points to form a closed curve and the curve is not already closed
    m_closedCurve = true;
    // Ensure the first and last control points are adjusted for smooth continuity
    UpdateTwoCtrlPoints(0);
    UpdateTwoCtrlPoints(anchorPoints.size()-1);
  }


  void drawAll() {
    // draw the whole line
    drawFromTo(0, anchorPoints.size() - (m_closedCurve ? 1 : 2));
  }
  
  void drawAllTranslate(PVector transVec) {
    // draw the whole line
    drawFromToTranslate(0, anchorPoints.size() - (m_closedCurve ? 1 : 2), transVec);
  }
  

  void drawFrom(int startIndex) {
    // draw from a certain point to the end
    drawFromTo(startIndex, anchorPoints.size() - (m_closedCurve ? 1 : 2));
  }
  
  void drawFromTo(int startIndex, int endIndex) {
    // draw some of the segments
    noFill();
    startIndex = max(startIndex, 0);
    endIndex = min(endIndex, anchorPoints.size() - (m_closedCurve ? 1 : 2));
    int nextIdx;
    int pointsCnt = anchorPoints.size();
    
    for (int i = startIndex; i <= endIndex; i++) {
      nextIdx = (i + 1) % pointsCnt; // in case it wraps after the last point (in a closed curve)
      bezier(anchorPoints.get(i).x, anchorPoints.get(i).y,
             ctrlPointsR.get(i).x, ctrlPointsR.get(i).y,
             ctrlPointsL.get(nextIdx).x, ctrlPointsL.get(nextIdx).y,
             anchorPoints.get(nextIdx).x, anchorPoints.get(nextIdx).y);
    }
  }

  void drawFromToTranslate(int startIndex, int endIndex, PVector transVec) {
    noFill();
    startIndex = max(startIndex, 0);
    endIndex = min(endIndex, anchorPoints.size() - (m_closedCurve ? 1 : 2));
    int nextIdx;
    int pointsCnt = anchorPoints.size();
    
    for (int i = startIndex; i <= endIndex; i++) {
      nextIdx = (i + 1) % pointsCnt; // in case it wraps after the last point (in a closed curve)
      bezier(anchorPoints.get(i).x + transVec.x, anchorPoints.get(i).y + transVec.y,
             ctrlPointsR.get(i).x + transVec.x, ctrlPointsR.get(i).y + transVec.y,
             ctrlPointsL.get(nextIdx).x + transVec.x, ctrlPointsL.get(nextIdx).y + transVec.y,
             anchorPoints.get(nextIdx).x + transVec.x, anchorPoints.get(nextIdx).y + transVec.y);
    }
  }
  
  void drawAnchorPoints() {
    int savedStrokeColor = g.strokeColor;
    stroke (0,0,255);
    for (PVector point : anchorPoints) {
      circle(point.x, point.y,5);
    }
    stroke(savedStrokeColor);
  }
 
  private void updateControlPoints() {
    // this is the smart part of the class 
    // the bezier segments are C1 continuous
    // after you have added one or several anchor-points, this procedure is called
    // to calculate the new control points that are still missing  and adjust some of the existing ones
    
    // here we find out, which points have no associated control points yet, i.e. where to start with the calculation
    // we don't want to start all over everytime some points are added
    // however, we do have to correct the control point that belongs to the last point, because now we know how
    // the line continues. 

    int startIdx = ctrlPointsL.size();               // how many control points are there currently 
    int anchorPointsCnt = anchorPoints.size();       // how many anchor points are there
    
    // add all the control points that are missing. Set them to 0 initially and update them later.
    for (int i = startIdx; i < anchorPointsCnt; i++) {
      ctrlPointsL.add(new PVector(0,0));
      ctrlPointsR.add(new PVector(0,0));
    }

    // if there is only one point, then we are done, since it does not need control points calculated
    if (anchorPointsCnt == 1) {
      return;
    }
    
    // now calculate or update all the affected control points, i.e. all the new controlPoints 
    // plus one to the left of the newly added anchor points (formerly the last control point)
    for (int i = max(startIdx - 1, 0); i < anchorPointsCnt; i++) {
      UpdateTwoCtrlPoints(i);
    }
  }
  
  private void UpdateTwoCtrlPoints(int i) {
    // Update both control points that belong to anchorPoint[i]
        
    // get the direction of the line connecting the point before and after this point.
    // if it is the first or last point in the line, it is just the direction of the following
    // or previous point respectively.
    // We always assume, that the line is not closed when points are added
    PVector tangent;
    
    
    int previousPointIdx;
    int nextPointIdx;
    
    if (m_closedCurve) {
      // a closed curve
      int pointsCnt = anchorPoints.size();
      previousPointIdx = ((i - 1) % pointsCnt + pointsCnt) % pointsCnt; // wrap around if necessary, ensure index is positive
      nextPointIdx = (i + 1) % pointsCnt;                               // wrap around if necessary
    } else {
      // an open line
      // handle the special cases of the first point and the last point in an open line
      previousPointIdx = max(0, i - 1);
      nextPointIdx = min(anchorPoints.size() - 1, i + 1);
    }
    
    // get the tangent line between the previous point and the next point to this point
    tangent = PVector.sub(anchorPoints.get(nextPointIdx), anchorPoints.get(previousPointIdx));
    // scale the line to length = 1
    tangent.normalize();
  
    // for the line to be C1 continues, the anchor-point and it's control points to the left
    // and right of it need to be co-linear (i.e. be on one straight line)
    // Now we calculate how far away form the point the two control points should be (while still on that imaginary straight line)
    
        
    float controlDistLeft = anchorPoints.get(i).dist(anchorPoints.get(previousPointIdx)) * m_swell;
    float controlDistRight = anchorPoints.get(i).dist(anchorPoints.get(nextPointIdx)) * m_swell;
    
    ctrlPointsL.set(i, PVector.sub(anchorPoints.get(i), tangent.copy().mult(controlDistLeft)));
    ctrlPointsR.set(i, PVector.add(anchorPoints.get(i), tangent.copy().mult(controlDistRight)));

  }
  
  public String formatNumber(float num) {
    // 3 digits after the comma should be plenty enough for a pen plotter
    return String.format(Locale.US, "%.3f", num);
  }
  
  String toSVGPath() {
    return toSVGPath(0,0,0);
  }
  
  String toSVGPath(float maxLength, int clipStart, int clipEnd) {
    // instead of using beginRecord and endRecord we use our own manual method to create an SVG-file.
    // That way, instead or recording a path for each bezier segement, we can make it one long
    // bezier path and ensure, that it will be drawn without a single pen lift.
    // If it is just a bunch of small paths, you never know how the plotter-SW will treat it.

    // Additionally we can split the long line up into several paths of given length (assuming 1 unit = 1 mm),
    // each path with a given max. length (or a bit over it)
    // This will make it much easier later to plot the line in several installments and re-fill 
    // the ink in the pen in between paths.
    
    // Additionally (added later); we provide a way to clip a few points at the start and the end.
    // Set both params to 0 if you want to write the full line.
        
    StringBuilder pathData = new StringBuilder();
    float totalLength = 0;
    int pointsCnt = anchorPoints.size();
    int nextIdx;
    
    if (pointsCnt < 2) return ""; // Return empty string if not enough points
    
    pathData.append("<path d=\"");

    // Start with the "move to" command for the first point
    PVector start = anchorPoints.get(0 + clipStart);
    pathData.append("M" + formatNumber(start.x) + " " + formatNumber(start.y));
  
    // Iterate through the list of points and append the Bezier curve commands
    for (int i = 0 + clipStart; i <= (pointsCnt - (m_closedCurve ? 1 : 2) - clipEnd); i++) {
      nextIdx = (i + 1) % pointsCnt; // in case it wraps
      PVector cp1 = ctrlPointsR.get(i);
      PVector cp2 = ctrlPointsL.get(nextIdx);
      PVector end = anchorPoints.get(nextIdx);
  
      pathData.append(" C" + formatNumber(cp1.x) + " " + formatNumber(cp1.y) + ", " +
                      formatNumber(cp2.x) + " " + formatNumber(cp2.y) + ", " +
                      formatNumber(end.x) + " " + formatNumber(end.y));
                      
      if (maxLength > 0) { // if the parameter is set to 0, we make only one path
        // calculate the path length so far
        totalLength += 
        bezierLength(anchorPoints.get(i), ctrlPointsR.get(i), ctrlPointsL.get(nextIdx), anchorPoints.get(nextIdx), 5);
        
        if (totalLength >= maxLength) {
          // end this path and start a new one
          pathData.append("\" fill=\"none\" stroke=\"black\" stroke-width=\"0.1\" />\n"); 
          // 
          pathData.append("<path d=\"");
          pathData.append("M" + formatNumber(end.x) + " " + formatNumber(end.y));
          totalLength = 0;
        }
      }
    }
  
    pathData.append("\" fill=\"none\" stroke=\"black\" stroke-width=\"0.1\" />\n"); 

    return pathData.toString();
  }
    
  float getLength(int samplesPerSegment) {
    // return the length of the whole line. Usually about 5 samplesPerSegment is precise enough
    return getLength(0, anchorPoints.size() - (m_closedCurve ? 1 : 2), samplesPerSegment);
  }
  
  float getLength(int startIndex, int endIndex, int samplesPerSegment) {
    // we approximate the total length by approximating the bezier curve with
    // a number of straight lines. The more samples we take, the more precise it gets (but also slower)
    float totalLength = 0.0;
    int pointsCnt = anchorPoints.size();
    
    startIndex = max(startIndex, 0);
    endIndex = min(endIndex, anchorPoints.size() - (m_closedCurve ? 1 : 2));
    for (int i = startIndex; i <= endIndex; i++) {
      totalLength += 
      bezierLength(anchorPoints.get(i), ctrlPointsR.get(i), ctrlPointsL.get((i+1) % pointsCnt), anchorPoints.get((i+1) % pointsCnt), samplesPerSegment);
    }
    return totalLength;
  }
  

  PVector getPointAtPct(float distPct) {
    // returns the point at a certain distance from the start of the line
    // the distance is given as a fraction of the total length and has to be between 0.0 and 1.0 (a percentage)
    
    // Ensure the percentage is between 0.0 and 1.0
    distPct = constrain(distPct, 0.0, 1.0);
  
    // Calculate the total length of the entire Bezier line
    float totalLength = getLength(100); // Increase samples for more accuracy
  
    // Calculate the target length along the curve based on the percentage
    float targetLength = distPct * totalLength;
  
    return getPointAtDist(targetLength);
  }


  PVector getPointAtDist(float dist) {
    // returns the point at a certain distance from the start of the line
    // the distance is given in mm and has to be between 0.0 amd the total length of the line

    // Variables to accumulate the length and track the previous point
    float accumulatedLength = 0.0;
    PVector prevPoint = anchorPoints.get(0);
    int pointsCnt = anchorPoints.size();
  
    // Iterate over each segment and calculate lengths more accurately
    for (int i = 0; i <= pointsCnt - (m_closedCurve ? 1 : 2); i++) {
      PVector p0 = anchorPoints.get(i);
      PVector p1 = ctrlPointsR.get(i);
      PVector p2 = ctrlPointsL.get((i + 1) % pointsCnt);
      PVector p3 = anchorPoints.get((i + 1) % pointsCnt);
  
      // Accumulate length along the segment, checking each small sub-segment
      int numSamples = 100; // More samples within each segment for precision
      for (int j = 1; j <= numSamples; j++) {
        float t = j / (float) numSamples;
        PVector currentPoint = bezierPoint(p0, p1, p2, p3, t);
        float segmentLength = prevPoint.dist(currentPoint);
  
        if (accumulatedLength + segmentLength >= dist) {
          // Interpolate within this small sub-segment
          float remainingLength = dist - accumulatedLength;
          float interpolationFactor = remainingLength / segmentLength;
          return PVector.lerp(prevPoint, currentPoint, interpolationFactor);
        }
  
        accumulatedLength += segmentLength;
        prevPoint = currentPoint;
      }
    }
  
    // If we reach here, return the last point on the curve (distPct was 1.0)
    return (m_closedCurve ? anchorPoints.get(0) : anchorPoints.get(pointsCnt - 1));    
  }

  PVector[] getEquidistantPointArr(int numPoints) {
    // returns an array of points that are distributed evenly over the whole length of the line
    // if the line is closed, it will return the same points as an open line, i.e. it will not return the first point twice
    
    // Calculate the total length of the bezier curve
    float totalLength = getLength(100);  // Or another appropriate number of samples per segment
    
    // Determine the distance between each point
    float distanceBetweenPoints = totalLength / (numPoints - (m_closedCurve ? 0 : 1));
  
    // Array to store the equidistant points
    PVector[] equidistantPoints = new PVector[numPoints];
    equidistantPoints[0] = anchorPoints.get(0);  // Start at the first point
  
    // Variables to keep track of progress along the curve
    float accumulatedLength = 0.0;
    PVector prevPoint = anchorPoints.get(0);
    int pointsCnt = anchorPoints.size();
    int currentPointIndex = 1;
  
    // Loop over each segment to find equidistant points
    for (int i = 0; i <= pointsCnt - (m_closedCurve ? 1 : 2); i++) {
      PVector p0 = anchorPoints.get(i);
      PVector p1 = ctrlPointsR.get(i);
      PVector p2 = ctrlPointsL.get((i + 1) % pointsCnt);
      PVector p3 = anchorPoints.get((i + 1) % pointsCnt);
  
      int numSamples = 100; // Use a constant number of samples within each segment for precision
      for (int j = 1; j <= numSamples; j++) {
        float t = j / (float) numSamples;
        PVector currentPoint = bezierPoint(p0, p1, p2, p3, t);
        float segmentLength = prevPoint.dist(currentPoint);
        
        // Check if the next equidistant point should be placed within this segment
        while (accumulatedLength + segmentLength >= currentPointIndex * distanceBetweenPoints && currentPointIndex < numPoints) {
          float remainingLength = currentPointIndex * distanceBetweenPoints - accumulatedLength;
          float interpolationFactor = remainingLength / segmentLength;
          equidistantPoints[currentPointIndex] = PVector.lerp(prevPoint, currentPoint, interpolationFactor);
          currentPointIndex++;
        }
  
        // Update accumulated length and move to the next segment
        accumulatedLength += segmentLength;
        prevPoint = currentPoint;
      }
    }
  
    // Ensure the last point is the end of the curve (unless it´s a closed curve, then we should never get here)
    if (currentPointIndex < numPoints && !m_closedCurve) {
      equidistantPoints[currentPointIndex] = anchorPoints.get(pointsCnt - 1);
    }
  
    return equidistantPoints;
  }
  
}  
  
// End of the class
// Helper function needed often in the class

float bezierLength(PVector p0, PVector p1, PVector p2, PVector p3, int numSamples) {
  // calculate the length of a single bezier segment
  float totalLength = 0;
  PVector prevPoint = p0;
  
  for (int i = 1; i <= numSamples; i++) {
    float t = i / (float)numSamples;
    PVector currentPoint = bezierPoint(p0, p1, p2, p3, t);
    totalLength += prevPoint.dist(currentPoint);
    prevPoint = currentPoint;
  }
  
  return totalLength;
}

PVector bezierPoint(PVector p0, PVector p1, PVector p2, PVector p3, float t) {
  // get the point somewhere on a bezier segment. 0 = start. 1 = end of line.
  float u = 1 - t;
  float tt = t*t;
  float uu = u*u;
  float uuu = uu * u;
  float ttt = tt * t;

  PVector result = p0.copy().mult(uuu);
  result.add(p1.copy().mult(3 * uu * t));
  result.add(p2.copy().mult(3 * u * tt));
  result.add(p3.copy().mult(ttt));

  return result;
}

// here come some more helper functions for writing to a svg file

void saveAsSVG(BezierLine line) {
  // create a new file containg the whole line in one path
  String svgContent = createSVGHeader(width, height);
  
  // add the bezier line path
  svgContent += line.toSVGPath();
  // add the footer
  svgContent += createSVGFooter();
  
  // write it to a new file in the current dorectory
  String timestamp = year() + "-" + month() + "-" + day() + "_" + hour() + "-" + minute() + "-" + second();
  saveSVGToFile(svgContent.toString(), sketchPath("output_" + timestamp + ".svg"));
  println(sketchPath("output_" + timestamp + ".svg") + " written");
}

String createSVGHeader(int width, int height) {
  // Create the header of an SVG-File
  StringBuilder svgContent = new StringBuilder();
  // 1. SVG Header
  svgContent.append("<?xml version=\"1.0\" standalone=\"no\"?>\n");
  svgContent.append("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
  svgContent.append("\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");

  // 2. Opening SVG tag (set the width and height to the canvas or desired dimensions)
  svgContent.append(String.format("<svg width=\"%d\" height=\"%d\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n", width, height));
  
  return svgContent.toString();
}

String createSVGHeader(int width, int height, String[] comments) {
  // Create the header of an SVG-File with some comments
  StringBuilder svgContent = new StringBuilder();
  // 1. SVG Header
  svgContent.append("<?xml version=\"1.0\" standalone=\"no\"?>\n");
  svgContent.append("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
  svgContent.append("\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");

  // 2. Comments
  if (comments != null && comments.length > 0) {
    svgContent.append("<!--\n");
    for (String comment : comments) {
        svgContent.append(comment).append("\n");
    }
    svgContent.append("-->\n");
  }

  // 3. Opening SVG tag (set the width and height to the canvas or desired dimensions)
  svgContent.append(String.format("<svg width=\"%d\" height=\"%d\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n", width, height));

  return svgContent.toString();
}


String createSVGFooter() {
  return ("</svg>");
}  

void saveSVGToFile(String svgContent, String filePath) {
  try {
    PrintWriter writer = new PrintWriter(filePath, "UTF-8");
    writer.print(svgContent);
    writer.close();
  } catch(Exception e) {
    e.printStackTrace();
  }
}


/*

// TESTING and USAGE EXAMPLE
// Every time you click left, a number of random points get added to the line
// Every time you click right, the whole line gets saved to a file and the line re-drawn
// Note: Because we only add a number of points while drawing and only draw the new segments, 
//       the connection between the larger segments is not continuous. We would have to delete the
//       last bezier curve on the screen and redraw it to draw it correctly, but that is slow.
//       Instead we draw it all new when saving the file.
//       The line inside the class is always correct, it's just that the last existing segment does not know how the line will later
//       continue and can therefore not end in the correct angle on the screen. As soon as you add points, this gets corrected.

final int N = 10;
BezierLine line;

void setup() {
  size(800, 800);
  background(255);
  line = new BezierLine();
}

void draw() {
  // Our animation loop will be mostly empty as we're just waiting for mouse input.
}

void mouseClicked() {

  if (mouseButton == LEFT) {
    // Create a list to store the new random points
    PVector[] newPoints = new PVector[N];
    
    for (int i = 0; i < N; i++) {
        newPoints[i] = new PVector(50 + random(width - 100), 50 + random(height - 100));
    }
    line.addPoints(newPoints);
    // only draw the new ones
    line.drawFromTo(line.anchorPoints.size() - N - 1, line.anchorPoints.size()-2);  
    println(line.getLength(5)/1000.0);

  } else {

    background(255);
    line.drawAll();
    // draw it anew and make sure all control points are calculated
    
    // save the long bezier line in an SVG in a single path
    // create the file header
    String svgContent = createSVGHeader(width, height);
    // add the bezier line path
    svgContent += line.toSVGPath();
    // add the footer
    svgContent += createSVGFooter();
    
    // write it a file in the current dorectory
    String timestamp = year() + "-" + month() + "-" + day() + "_" + hour() + "-" + minute() + "-" + second();
    saveSVGToFile(svgContent.toString(), sketchPath("output_" + timestamp + ".svg"));
    println(sketchPath("output_" + timestamp + ".svg") + " written");
    
    println(line.getLength(5)/1000.0);

  }
}

*/

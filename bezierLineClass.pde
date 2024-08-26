// Author: matt@void.ch
// This class represents a C1 continuous, smooth Bezier curve through all the points that you give it
// The data contains the anchor-points (given by you) and the control-points (calculated by the class)

import java.util.Locale;

class BezierLine {
  ArrayList<PVector> points;
  ArrayList<PVector> controlPoints;
  
  BezierLine() {
    points = new ArrayList<PVector>();
    controlPoints = new ArrayList<PVector>();
  }

  void addPoint(PVector point) {
    // add a single anchor point 
    points.add(point);
    calculateControlPoints();  // Will implement this in a moment
  }
  
  void movePoint(int index, PVector moveVector) {
    // move a single anchor point by a vector 
    points.get(index).x += moveVector.x;
    points.get(index).y += moveVector.y;
    // recalc everything (it's a bit overkill: todo: only re-calc the controlpoints necessary
    controlPoints.clear();
    calculateControlPoints(); 
  }

  void addPoints(PVector[] newPoints) {
    // add an array of points 
    for (PVector p : newPoints) {
      points.add(p);
    }
    calculateControlPoints();  // Adjust to only calculate necessary control points
  }
  
  void addPoints(ArrayList<PVector> newPointsList) {
    // add a list of points 
    for (PVector p : newPointsList) {
      points.add(p);
    }
    calculateControlPoints();  // Adjust to only calculate necessary control points
  }
  
  void drawAll() {
    // draw the whole line
    drawFromTo(0, points.size()-2);
  }
  
  void drawAllTranslate(PVector transVec) {
    // draw the whole line
    drawFromToTranslate(0, points.size()-2, transVec);
  }
  

  void drawFrom(int startIndex) {
    // draw from a certain point to the end
    drawFromTo(startIndex, points.size()-2);
  }
  
  void drawFromTo(int startIndex, int endIndex) {
    // draw some of the segments
    noFill();
    startIndex = max(startIndex, 0);
    endIndex = min(endIndex, points.size() - 2);
    for (int i = startIndex; i <= endIndex; i++) {
        /*
        // temp hack to try something
        if ((i % 3) == 0) {
          stroke(255,0,0);
        }
        if ((i % 3) == 1) {
          stroke(0,255,0);
        }
        */
        bezier(points.get(i).x, points.get(i).y,
               controlPoints.get(2*i).x, controlPoints.get(2*i).y,
               controlPoints.get(2*i+1).x, controlPoints.get(2*i+1).y,
               points.get(i+1).x, points.get(i+1).y);
        // cont. hack
        // stroke(0,0,0);
        
    }
  }

  void drawFromToTranslate(int startIndex, int endIndex, PVector transVec) {
    // draw some of the segments
    noFill();
    startIndex = max(startIndex, 0);
    endIndex = min(endIndex, points.size() - 2);
    for (int i = startIndex; i <= endIndex; i++) {
        bezier(points.get(i).x+transVec.x, points.get(i).y+transVec.y,
               controlPoints.get(2*i).x+transVec.x, controlPoints.get(2*i).y+transVec.y,
               controlPoints.get(2*i+1).x+transVec.x, controlPoints.get(2*i+1).y+transVec.y,
               points.get(i+1).x+transVec.x, points.get(i+1).y+transVec.y);
    }
  }
 
 
  private void calculateControlPoints() {
    // this is the smart part of the class 
    // the bezier segments are C1 continuous
    // after you have added one or several anchor-points, this procedure is called
    // to calculate the new control points that are still missing
    
    if (points.size() < 2) return;
  
    // here we find out, which points have no associated control points yet, i.e. where to start with the calculation
    // we don't want to start all over everytime some points are added
    // however, we do have to correct the control point that belongs to the last point, because now we know how
    // the line continues. This makes the code a little messy

    int startIdx = controlPoints.size() / 2;
  
    for (int i = startIdx; i < points.size(); i++) {
        // get the direction of the line connecting the point before and after this point.
        // if it is the first or last point in the line, it is just the direction of the following
        // or previous point respectively
        PVector tangent;
        if (i == 0) {
            // first point
            tangent = PVector.sub(points.get(1), points.get(0));
        } else if (i == points.size() - 1) {
            // last point
            tangent = PVector.sub(points.get(i), points.get(i-1));
        } else {
            // all the points in between
            tangent = PVector.sub(points.get(i+1), points.get(i-1));
        }
        // scale the line to length = 1
        tangent.normalize();
  
        // for the line to be C1 continues, the anchot-point and it's control points to the left
        // and right of it need to be co-linear (i.e. be on one straight line)
        // Now we calculate how far away form the point the two control points should be (while still on that imaginary straight line)
        // todo: we could make this hard-coded factor of 0.3 a parameter and therefore change the "curviness" of the line
        
        float controlDistLeft, controlDistRight;
  
        if (i == 0) {
            // first point in the whole line; there is only a control point to the right of it
            controlDistRight = points.get(i).dist(points.get(i+1)) * 0.3f;
            controlPoints.add(PVector.add(points.get(i), tangent.copy().mult(controlDistRight)));
        } else if (i == points.size() - 1) {
            // last point in the whole line; there is only a control point to the left of it
            controlDistLeft = points.get(i).dist(points.get(i-1)) * 0.3f;
            controlPoints.add(PVector.sub(points.get(i), tangent.copy().mult(controlDistLeft)));
        } else {
            // all the other points in between
            controlDistLeft = points.get(i).dist(points.get(i-1)) * 0.3f;
            controlDistRight = points.get(i).dist(points.get(i+1)) * 0.3f;
  
            // If it's the starting index of the new points, use `set` for correcting the already existing left control point
            // otherwise add the new control-points to the left
            if (i == startIdx) {
                controlPoints.set(2*(i-1)+1, PVector.sub(points.get(i), tangent.copy().mult(controlDistLeft)));
            } else {
                controlPoints.add(PVector.sub(points.get(i), tangent.copy().mult(controlDistLeft)));
            }
            // now add the right control-point in both cases
            controlPoints.add(PVector.add(points.get(i), tangent.copy().mult(controlDistRight)));

        }
    }
    // println(points.size(), controlPoints.size());
  }
  
  public String formatNumber(float num) {
      return String.format(Locale.US, "%.4f", num).replaceAll("0*$", "").replaceAll("\\.$", "");
  }

  String toSVGPath() {
    // instead of using beginRecord and endRecord we use our own manual method to create an SVG-file.
    // That way, instead or recording a path for each bezier segement, we can make it one long
    // bezier path and insure, that it will be drawn without a single pen lift.
    // if it is a bunch of paths, you never know how the plotte-SW will treat it
    // also: Later, is we wish, we can expand this saving-function to create a new layer every X meter of
    // line-length, therefore enabling us to print layer by layer and refill the ink or change the pen
    // before they run out. Many pens only manage about 150 to 250 meters before they run out.
    // also, we can control the precision. 2 digits after the comma should be plenty enough for 
    // a pen plotter and keeps the file smaller
    
    StringBuilder pathData = new StringBuilder();
    
    if (points.size() < 2) return ""; // Return empty string if not enough points
    
    pathData.append("<path d=\"");

    // Start with the "move to" command for the first point
    PVector start = points.get(0);
    pathData.append("M" + formatNumber(start.x) + " " + formatNumber(start.y));
  
    // Iterate through the list of points and append the Bezier curve commands
    for (int i = 0; i < points.size() - 1; i++) {
      PVector cp1 = controlPoints.get(2 * i);
      PVector cp2 = controlPoints.get(2 * i + 1);
      PVector end = points.get(i + 1);
  
      pathData.append(" C" + formatNumber(cp1.x) + " " + formatNumber(cp1.y) + ", " +
                      formatNumber(cp2.x) + " " + formatNumber(cp2.y) + ", " +
                      formatNumber(end.x) + " " + formatNumber(end.y));
    }
  
    pathData.append("\" fill=\"none\" stroke=\"black\" />\n");

    return pathData.toString();
  }
  
  String toSVGPathsMaxLengthPerPath(float maxLength, int clipStart, int clipEnd) {
    // as above, but additionally we split the long line up several in paths,
    // each path with a given max. length (or a bit over it)
    // This will make it much easier later, to print the line in several installments and re-fill 
    // the ink in the pen in between paths
    // Additionally we provide a way to clip a few points at the start and the end
    // set both params to 0 if you want to write the full line

        
    StringBuilder pathData = new StringBuilder();
    float totalLength = 0;
    
    if (points.size() < 2) return ""; // Return empty string if not enough points
    
    pathData.append("<path d=\"");

    // Start with the "move to" command for the first point
    PVector start = points.get(0 + clipStart);
    pathData.append("M" + formatNumber(start.x) + " " + formatNumber(start.y));
  
    // Iterate through the list of points and append the Bezier curve commands
    for (int i = 0 + clipStart; i < (points.size() - 1 - clipEnd); i++) {
      PVector cp1 = controlPoints.get(2 * i);
      PVector cp2 = controlPoints.get(2 * i + 1);
      PVector end = points.get(i + 1);
  
      pathData.append(" C" + formatNumber(cp1.x) + " " + formatNumber(cp1.y) + ", " +
                      formatNumber(cp2.x) + " " + formatNumber(cp2.y) + ", " +
                      formatNumber(end.x) + " " + formatNumber(end.y));
                      
      // calculate the path length so far
      totalLength += 
      bezierLength(points.get(i), controlPoints.get(2*i), controlPoints.get(2*i+1), points.get(i+1), 5);
      
      if (totalLength >= maxLength) {
        // end this path and start a new one
        pathData.append("\" fill=\"none\" stroke=\"black\" stroke-width=\"0.1\" />\n"); 
        // 
        pathData.append("<path d=\"");
        pathData.append("M" + formatNumber(end.x) + " " + formatNumber(end.y));
        totalLength = 0;
      }
    }
  
    pathData.append("\" fill=\"none\" stroke=\"black\" stroke-width=\"0.1\" />\n"); 

    return pathData.toString();
  }
    
  
  float getLenght(int startIndex, int endIndex, int samplesPerSegment) {
    // we approximate the total length by approximating the bezier curve with
    // a number of straight lines. The more samples we take, the more precise it gets (but also slower).
    float totalLength = 0.0;
    startIndex = max(startIndex, 0);
    endIndex = min(endIndex, points.size() - 2);
    for (int i = startIndex; i <= endIndex; i++) {
      totalLength += 
      bezierLength(points.get(i), controlPoints.get(2*i), controlPoints.get(2*i+1), points.get(i+1), samplesPerSegment);
    }
    return totalLength;
  }
  
  PVector getPointAtPct(float distPct) {
    // Ensure the percentage is between 0.0 and 1.0
    distPct = constrain(distPct, 0.0, 1.0);
  
    // Calculate the total length of the entire Bezier line
    float totalLength = getLenght(0, points.size() - 2, 100); // Increase samples for more accuracy
  
    // Calculate the target length along the curve based on the percentage
    float targetLength = distPct * totalLength;
  
    // Variables to accumulate the length and track the previous point
    float accumulatedLength = 0.0;
    PVector prevPoint = points.get(0);
  
    // Iterate over each segment and calculate lengths more accurately
    for (int i = 0; i < points.size() - 1; i++) {
      PVector p0 = points.get(i);
      PVector p1 = controlPoints.get(2 * i);
      PVector p2 = controlPoints.get(2 * i + 1);
      PVector p3 = points.get(i + 1);
  
      // Accumulate length along the segment, checking each small sub-segment
      int numSamples = 100; // More samples within each segment for precision
      for (int j = 1; j <= numSamples; j++) {
        float t = j / (float) numSamples;
        PVector currentPoint = bezierPoint(p0, p1, p2, p3, t);
        float segmentLength = prevPoint.dist(currentPoint);
  
        if (accumulatedLength + segmentLength >= targetLength) {
          // Interpolate within this small sub-segment
          float remainingLength = targetLength - accumulatedLength;
          float interpolationFactor = remainingLength / segmentLength;
          return PVector.lerp(prevPoint, currentPoint, interpolationFactor);
        }
  
        accumulatedLength += segmentLength;
        prevPoint = currentPoint;
      }
    }
  
    // If we reach here, return the last point on the curve (distPct was 1.0)
    return points.get(points.size() - 1);
  }

  PVector getPointAtDist(float dist) {
    // Calculate the total length of the entire Bezier line
    float totalLength = getLenght(0, points.size() - 2, 100); // 100 samples for more accuracy
  
    // Ensure the distance is within the valid range
    dist = constrain(dist, 0.0, totalLength);
  
    // Calculate the percentage corresponding to the given distance
    float distPct = dist / totalLength;
  
    // Call the existing function to get the point at the calculated percentage
    return getPointAtPct(distPct);
  }  
}

// End of the class
// here come some more helper functions

String createSVGHeader(int width, int height) {
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
  StringBuilder svgContent = new StringBuilder();
  svgContent.append("</svg>");
  
  return svgContent.toString();
}  

void saveSVGToFile(String svgContent, String filePath) {
  try {
    PrintWriter writer = new PrintWriter(filePath, "UTF-8");
    writer.print(svgContent);
    writer.close();
    println ("written");
  } catch(Exception e) {
    e.printStackTrace();
  }
}

// maybe use later to calculate the total length of a line
// and to be able to chop up the SVG-output every X meters
// and maybe pack these bits into different inkscape-layers, so that 
// we can check and refill the ink every X meters

float bezierLength(PVector p0, PVector p1, PVector p2, PVector p3, int numSamples) {
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

/*

// TESTING
// Every time you click left a number of random points get added to the line
// Every time you click right, the whole line gets saved to a file a the line re-drawn
// Note: Because we only add a number of points while drawing and only draw the new segmenst, 
//       The connection between the larger segments is not continuous. We would have to delete the
//       last bezier curve on the screen and redraw it, but that is messy
//       instead we draw it all new when saving the file.
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
    line.drawFromTo(line.points.size() - N - 1, line.points.size()-2);  
    println(line.getLenght(0, line.points.size(), 5)/1000.0);

} else {
    // save the line in an SVG
    background(255);
    line.drawAll();
    // println(line.toSVGPath());
    
    String svgContent = createSVGHeader(width, height);
    svgContent += line.toSVGPath();
    svgContent += createSVGFooter();
    
    // println(svgContent);
    saveSVGToFile(svgContent, sketchPath("output.svg"));
    
    println(line.getLenght(0, line.points.size(), 5)/1000.0);

  }
}

*/

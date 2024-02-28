package info.quantlab.easyplot;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.*;
import java.util.function.Function;

public class PlotRunner {

    public static void printHelp() {
        String[] commands = new String[] {
                "saveasjava",
                "saveassvg", "saveaspng", "saveasjpg", "saveaspdf",
                "settitle",
                "setxaxislabel",
                "setyaxislabel",
                "setislegendvisible",
                "changeplotcolor",
                "changeplotname",
                "removeplot",
                "removemarker",
                "setnextmarker",
                "setplotstroke",
                "help"
        };
        for(String command: commands) {
            System.out.println(methodNameSignature(command));
        }
    }

    public static String methodNameSignature(String command) {
        String commandNoLowerCase;
        command = command.toLowerCase();
        return switch (command) {
            case "saveasjava" -> "SaveAsJava has 0 arguments.\n" +
                    "\tDescription: Saves a runnable Java file for the fast recreation of the plot. Cannot recover PlotStyle nor Size of the plot.";
            case "saveassvg", "saveaspng", "saveasjpg", "saveasjpeg", "saveaspdf" -> {
                commandNoLowerCase = command.substring(0, 1).toUpperCase() + command.substring(1, 4) + command.substring(4, 5).toUpperCase() + command.substring(5, 6) +
                        command.substring(6).toUpperCase();
                yield commandNoLowerCase + "has 3 Arguments:\n" +
                        "\t-file, -f, -path, -p or 0 as the filename/-path (String)\n" +
                        "\t-width, -w or 1 as the width of the picture (int)\n" +
                        "\t-height, -h or 2 as the height of the picture (int)\n" +
                        "\tDescription: Saves the plot as picture of the specified format.";
            }
            case "settitle" -> """
                    SetTitle has 1 argument.
                    \t-title, -t or 0 as the new title (String)
                    \tDescription: Changes the title of the plot.""";
            case "setxaxislabel" -> """
                    SetXAxisLabel has 1 argument.
                    \t-xAxisLabel, -l or 0 as the new xAxisLabel (String)
                    \tDescription: Changes the xAxisLabel of the plot.""";
            case "setyaxislabel" -> """
                    SetYAxisLabel has 1 argument.
                    \t-yAxisLabel, -l or 0 as the new yAxisLabel (String)
                    \tDescription: Changes the yAxisLabel of the plot.""";
            case "setislegendvisible" -> """
                    SetIsLegendvisible has 1 argument.
                    \t-isLegendvisible, -v or 0 as the new visibility of the legend
                    \tDescription: Changes visibility of the legend of the plot.""";
            case "changeplotcolor" -> """
                    ChangePlotColor has either 4 or 2 arguments.
                    \t-r, -red or 0 for red channel (int)
                    \t-g, -green or 1 for green channel (int)
                    \t-b, -blue or 2 for blue channel (int)
                    \t-i, -index or 3 for index of plotable to apply the color to (int)
                    "or one String and one Integer at:
                    \t-color or -colorName (String)
                    \t-i, -index or 1 for index of plotable to apply the color to (int)
                    \tDescription: Changes the plotcolor for the specified index. Specifying an invalid index or out of bounds color will throw an exception!""";
            case "changeplotname" -> """
                    ChangePlotName has 2 arguments.
                    \t-name or 0 as the new name (String)
                    \t-index or 1 as the index to change (int)
                    \tDescription: Changes the name of the plot corresponding to the index. Note this is the name of a single plotable in the plot.""";
            case "removeplot" -> """
                    RemovePlot has 1 argument.
                    \t-index, -i or 0 as the index to remove (int)
                    \tDescription: Removes the specified plot from the plot window. This is not reversible!""";
            case "removemarker" -> """
                    RemoveMarker has 1 argument.
                    \t-index, -i or 0 as the index to remove the markers from (int)
                    \tDescription: Removes the markers from the specified plot.""";
            case "setnextmarker" -> """
                    SetNextMarker has 1 argument.
                    \t-index, -i or 0 as the index to remove the markers from (int)
                    \tDescription: Sets the markers for the specified plot to be the next default markers.""";
            case "setplotstroke" -> """
                    SetPlotStroke has 2 arguments.
                    \t-index, -i or 0 as the index to remove the markers from (int)
                    \t-visibility, -v or 1 as the visibility of the stroke (boolean)
                    \tDescription: Sets the visibility of the stroke for the specified plot.""";
            case "help" -> """
                    Help has 0 arguments.
                    \tDescription: Prints helpful informations.""";
            default -> "Unknown Command!";
        };
    }

    public static void runConfigLoop(final EasyPlot2D[] plots) {
        Scanner scanner = new Scanner(System.in);
        String input;

        // Infinite loop

        while (true) {
            System.out.print("Enter a command (type 'exit' to quit): ");
            input = scanner.nextLine();

            if (input.equalsIgnoreCase("exit")) {
                System.out.println("Exiting the program...");
                break; // Exit the loop
            } else {
                String[] parts = input.split("\\s+");
                String command = parts[0].toLowerCase();

                Map<String, String> args = parseArguments(parts);

                if(args.containsKey("error")) {
                    System.out.println("Parsing the arguments yields an error! Error message: " + args.get("error"));
                    continue;
                }

                if(command.equals("help")) {
                    printHelp();
                    continue;
                }

                Function<Integer, Object> method = null;
                switch (command) {
                    case "saveasjava":
                        method = (index) -> plots[index].writeRecreationJavaFile();
                        break;
                    case "saveassvg":
                    case "saveaspng":
                    case "saveasjpg":
                    case "saveasjpeg":
                    case "saveaspdf":
                    {
                        String sFile = args.containsKey("0") ? args.get("0") :
                                (args.containsKey("f") ? args.get("f") : args.getOrDefault("file", null));
                        sFile = sFile != null ? sFile :
                                (args.containsKey("p") ? args.get("p") : args.getOrDefault("path", null));
                        String sWidth = args.containsKey("1")? args.get("1") :
                                (args.containsKey("w")? args.get("w") : args.getOrDefault("width", "1200"));
                        String sHeight = args.containsKey("2")? args.get("2") :
                                (args.containsKey("h")? args.get("h") : args.getOrDefault("height", "1200"));
                        try {
                            if(sFile == null || sFile.isEmpty()) {
                                throw new IllegalArgumentException();
                            }
                            final File file = new File(sFile);
                            final int width = Integer.parseInt(sWidth);
                            final int height = Integer.parseInt(sHeight);
                            switch(command) {
                                case "saveassvg":
                                    method = (index) -> {
                                        try {
                                            return plots[index].saveAsSVG(file, width, height);
                                        } catch (IOException e) {
                                            System.out.println("Saving the plot was unsuccessfull. Most likely the file is currently blocked or was given an illegal file name!");
                                            return new Object();
                                        }
                                    };
                                    break;
                                case "saveaspng":
                                    method = (index) -> {
                                        try {
                                            return plots[index].saveAsPNG(file, width, height);
                                        } catch (IOException e) {
                                            System.out.println("Saving the plot was unsuccessfull. Most likely the file is currently blocked or was given an illegal file name!");
                                            return new Object();
                                        }
                                    };
                                    break;
                                case "saveasjpg":
                                case "saveasjpeg":
                                    method = (index) -> {
                                        try {
                                            return plots[index].saveAsJPG(file, width, height);
                                        } catch (IOException e) {
                                            System.out.println("Saving the plot was unsuccessfull. Most likely the file is currently blocked or was given an illegal file name!");
                                            return new Object();
                                        }
                                    };
                                    break;
                                case "saveaspdf":
                                    method = (index) -> {
                                        try {
                                            return plots[index].saveAsPDF(file, width, height);
                                        } catch (IOException e) {
                                            System.out.println("Saving the plot was unsuccessfull. Most likely the file is currently blocked or was given an illegal file name!");
                                            return new Object();
                                        }
                                    };
                                    break;

                            }
                        }
                        catch (Exception ex){
                            System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                        }
                        break;
                    }
                    case "settitle":
                    {
                        final String title = args.containsKey("0") ? args.get("0") :
                                (args.containsKey("t") ? args.get("t") : args.getOrDefault("title", ""));
                        method = (index) -> plots[index].setTitle(title);
                        break;
                    }

                    case "setxaxislabel":
                    {
                        final String xAxisLabel = args.containsKey("0") ? args.get("0") :
                                (args.containsKey("l") ? args.get("l") : args.getOrDefault("xaxislabel", ""));
                        method = (index) -> plots[index].setXAxisLabel(xAxisLabel);
                        break;
                    }
                    case "setyaxislabel":
                    {
                        final String yAxisLabel = args.containsKey("0") ? args.get("0") :
                                (args.containsKey("l") ? args.get("l") : args.getOrDefault("yaxislabel", ""));
                        method = (index) -> plots[index].setYAxisLabel(yAxisLabel);
                        break;
                    }
                    case "setislegendvisible":
                    {
                        boolean isLegendVisible = false;
                        String visibilityArg = args.containsKey("0") ? args.get("0") :
                            (args.containsKey("v") ? args.get("v") : args.getOrDefault("islegendvisible", "false"));
                        try {
                            isLegendVisible = Boolean.parseBoolean(visibilityArg);
                        } catch (Exception ex){
                            try {
                                int intBool = Integer.parseInt(visibilityArg);
                                isLegendVisible = intBool != 0;
                            } catch (Exception ex2) {
                                System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                            }
                        }
                        final boolean legendArg = isLegendVisible;
                        method = (index) -> plots[index].setIsLegendVisible(legendArg);
                        break;
                    }
                    case "changeplotcolor":
                    {
                        try {
                            String red = args.containsKey("color") ? args.get("color") : args.getOrDefault("colorname", null);
                            String pIndex = null;
                            Color col = null;
                            if(red != null) {
                                pIndex = args.containsKey("1")? args.get("1") :
                                    (args.containsKey("i")? args.get("i") : args.getOrDefault("index", "0"));
                                Field field = Color.class.getField(red.toLowerCase());
                                col = (Color) field.get(null);
                            }
                            else {
                                String green;
                                String blue;
                                red = args.containsKey("0")? args.get("0") :
                                    (args.containsKey("r")? args.get("r") : args.getOrDefault("red", "0"));
                                green = args.containsKey("1")? args.get("1") :
                                    (args.containsKey("g")? args.get("g") : args.getOrDefault("green", "0"));
                                blue = args.containsKey("2")? args.get("2") :
                                    (args.containsKey("b")? args.get("b") : args.getOrDefault("blue", "0"));
                                pIndex = args.containsKey("3")? args.get("3") :
                                    (args.containsKey("i")? args.get("i") : args.getOrDefault("index", "0"));
                                int r = Integer.parseInt(red);
                                int g = Integer.parseInt(green);
                                int b = Integer.parseInt(blue);
                                col = new Color(r, g, b);
                            }

                            final int pInd = Integer.parseInt(pIndex);
                            final Color newColor = col;
                            method = (index) -> plots[index].changePlotColor(pInd, newColor);
                        } catch (Exception ex) {
                            System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                            break;
                        }
                    }
                    case "changeplotname":
                    {
                        final String name = args.containsKey("0") ? args.get("0") :
                                (args.containsKey("n") ? args.get("n") : args.getOrDefault("name", ""));
                        String pIndex = args.containsKey("1")? args.get("1") :
                                (args.containsKey("i")? args.get("i") : args.getOrDefault("index", "0"));
                        try {
                            final int plotIndex = Integer.parseInt(pIndex);
                            method = (index) -> plots[index].changePlotName(plotIndex, name);
                        } catch (Exception ex){
                            System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                        }
                        break;
                    }
                    case "removeplot":
                    {
                        String pIndex = args.containsKey("0")? args.get("0") :
                                (args.containsKey("i")? args.get("i") : args.getOrDefault("index", "-1"));
                        try {
                            final int plotIndex = Integer.parseInt(pIndex);
                            if(plotIndex < 0) {
                                throw new IllegalArgumentException();
                            }
                            method = (index) -> plots[index].removePlot(plotIndex);
                        } catch (Exception ex){
                            System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                        }
                        break;
                    }
                    case "removemarker":
                    {
                        String pIndex = args.containsKey("0")? args.get("0") :
                                (args.containsKey("i")? args.get("i") : args.getOrDefault("index", "0"));
                        try {
                            final int plotIndex = Integer.parseInt(pIndex);
                            if(plotIndex < 0) {
                                throw new IllegalArgumentException();
                            }
                            method = (index) -> plots[index].changePlotMarker(plotIndex, null);
                        } catch (Exception ex){
                            System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                        }
                        break;
                    }
                    case "setnextmarker":
                    {
                        String pIndex = args.containsKey("0")? args.get("0") :
                                (args.containsKey("i")? args.get("i") : args.getOrDefault("index", "0"));
                        try {
                            final int plotIndex = Integer.parseInt(pIndex);
                            if(plotIndex < 0) {
                                throw new IllegalArgumentException();
                            }
                            method = (index) -> plots[index].changePlotMarker(plotIndex, plots[index].skipShape());
                        } catch (Exception ex){
                            System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                        }
                        break;
                    }
                    case "setplotstroke":
                    {
                        String sIndex = args.containsKey("0")? args.get("0") :
                                (args.containsKey("i")? args.get("i") : args.getOrDefault("index", "-1"));
                        String sVis = args.containsKey("1")? args.get("1") :
                                (args.containsKey("v")? args.get("v") : args.getOrDefault("visible", "0"));
                        try {
                            final int plotIndex = Integer.parseInt(sIndex);
                            if(plotIndex < 0) {
                                throw new IllegalArgumentException();
                            }
                            boolean vis = false;
                            try {
                                vis = Boolean.parseBoolean(sVis);
                            } catch (Exception ex2) {
                                try {
                                    int intBool = Integer.parseInt(sVis);
                                    vis = intBool != 0;
                                } catch (Exception ex3) {
                                    throw new IllegalArgumentException();
                                }
                            }
                            final boolean visF = vis;
                            method = (index) -> plots[index].changePlotStroke(plotIndex, visF? new BasicStroke(2.0f) : null);
                        } catch (Exception ex){
                            System.out.println("Parsing the arguments yields an error!\n" + methodNameSignature(command));
                        }
                        break;
                    }
                    // Add more cases for other methods as needed
                    default:
                        break;

                }


                if (method == null) {
                    System.out.println("Method not found or invalid arguments! Input: " + input);
                    System.out.println("Type \"help\" for help");
                    continue;
                }
                int plotIndex = -1;
                try {
                    plotIndex = Integer.parseInt(args.getOrDefault("plotindex", "-1"));
                } catch(Exception ex) {
                    if(args.getOrDefault("plotindex", "-1").equalsIgnoreCase("a")) {
                        plotIndex = Integer.MAX_VALUE;
                    }
                }


                if(plots.length == 1) {
                    plotIndex = 0;
                }
                else if(plotIndex < 0) {
                    System.out.println("Which plot window should this method be applied to? Enter \"a\" for all plots.");
                    input = scanner.nextLine();
                    if(input.equalsIgnoreCase("a")) {
                        plotIndex = Integer.MAX_VALUE;
                    } else {
                        try {
                            plotIndex = Integer.parseInt(input);
                            if(plotIndex < 0) {
                                throw new IllegalArgumentException();
                            }
                        } catch(Exception ex) {
                            System.out.println("Did not work! Try another command!");
                            continue;
                        }
                    }
                }
                if(plotIndex == Integer.MAX_VALUE) {
                    for (int pIndex = 0; pIndex < plots.length; pIndex++) {
                        method.apply(pIndex);
                    }
                } else {
                    method.apply(plotIndex);
                }
            }


        }
        scanner.close(); // Close the scanner when done
    }

    public static Map<String, String> parseArguments(String[] args) {
        Map<String, String> argMap = new HashMap<>();
        ArrayList<String> realArgs = new ArrayList<>();
        ArrayList<String> argNames = new ArrayList<>();
        // Parse for quotes
        {
            // args[0] should always be the command
            int j=1;
            while (j < args.length) {
                String realArg = args[j];
                String argName = null;
                if (args[j].startsWith("\"") && args[j].endsWith("\"")) {
                    realArg = args[j].substring(1,args[j].length()-1);
                }
                else if(args[j].startsWith("\"")) {
                    StringBuilder stringBuilder = new StringBuilder(args[j].substring(1));
                    j++;
                    while(j < args.length) {
                        if(args[j].endsWith("\"")) {
                            stringBuilder.append(args[j], 0, args[j].length() - 1);
                            break;
                        }
                        stringBuilder.append(args[j]);
                        j++;
                    }
                    if(j == args.length) {
                        argMap.put("error", "Could not parse quotes!");
                        return argMap;
                    }
                    realArg = stringBuilder.toString();
                }
                else if (args[j].startsWith("-")) {
                    // Argument name found
                    argName = args[j].toLowerCase().substring(1);
                }
                if(argName == null) {
                    realArgs.add(realArg);
                    if(realArgs.size() > argNames.size()) {
                        argName = "" + argNames.size();
                        argNames.add(argName);
                        if(realArgs.size() > argNames.size()) {
                            argMap.put("error", "My Parser did not work correct! You should save the plots " +
                                    "and rerun with debug mode on");
                            return argMap;
                        }
                    }
                } else {
                    if(realArgs.size() < argNames.size()) {
                        argMap.put("error", "There are argument Names without Arguments!");
                        return argMap;
                    }
                    argNames.add(argName);
                }
                j++;
            }
            if(realArgs.size() != argNames.size()) {
                argMap.put("error", "My Parser did not work correct! You should save the plots " +
                        "and rerun with debug mode on");
                return argMap;
            }
        }
        for (int i = 0; i < realArgs.size(); i++) {
            argMap.put(argNames.get(i), realArgs.get(i));
        }
        return argMap;
    }
}

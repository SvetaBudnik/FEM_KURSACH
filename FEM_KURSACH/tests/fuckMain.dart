import 'dart:io';

class Node {
  double x = 0;
  double y = 0;

  Node({this.x = 0.0, this.y = 0.0});

  @override
  String toString() {
    return "(x: ${x}, y: ${y})";
  }

  String toFileString() {
    return "${x.toStringAsFixed(4).padLeft(8)} ${y.toStringAsFixed(4).padLeft(8)}";
  }
}

class FiniteRect {
  List<int> inds = [];

  FiniteRect({required this.inds});

  @override
  String toString() {
    return inds.toString();
  }

  String toFileString() {
    var out = "";
    for (var ind in inds) {
      out += "${ind.toString().padLeft(4)} ";
    }
    out += "  0";
    return out;
  }
}

class FiniteTreug {
  List<int> inds = [];

  FiniteTreug({required this.inds});

  @override
  String toString() {
    return inds.toString();
  }

  String toFileString() {
    var out = "";
    for (var ind in inds) {
      out += "${ind.toString().padLeft(4)} ";
    }
    out += "   0";
    return out;
  }
}

(int, List<Node>) CreateNodes(double beg, double end, double step) {
  var coords = <double>[];
  for (var i = beg; i <= end; i += step) {
    coords += [i];
  }

  var nodes = <Node>[];
  for (var i = 0; i < coords.length; i++) {
    for (var j = 0; j < coords.length; j++) {
      var node = Node(
        x: coords[j],
        y: coords[i],
      );
      nodes += [node];
    }
  }

  return (coords.length, nodes);
}

List<FiniteRect> GenerateFinites(int nodesInRow) {
  var finites = <FiniteRect>[];
  var finitesCount = (nodesInRow - 1) ~/ 2;

  for (int i = 0; i < finitesCount; i++) {
    for (int j = 0; j < finitesCount; j++) {
      var inds = <int>[];
      var firstInd = j * 2 + i * nodesInRow * 2;
      for (int k = 0; k < 3; k++) {
        for (int m = 0; m < 3; m++) {
          inds += [firstInd + m + k * nodesInRow];
        }
      }

      var finite = FiniteRect(inds: inds);
      finites += [finite];
    }
  }

  return finites;
}

List<FiniteTreug> GenerateTreugs(List<FiniteRect> rects) {
  var finites = <FiniteTreug>[];
  for (var rect in rects) {
    var i = rect.inds;
    var inds = [i[0], i[6], i[2], i[3], i[4], i[1]];
    finites.add(FiniteTreug(inds: inds));
    inds = [i[2], i[6], i[8], i[4], i[7], i[5]];
    finites.add(FiniteTreug(inds: inds));
  }

  return finites;
}

List<int> GenerateS1Nodes(int nodesInRow) {
  var s1Nodes = <int>[];
  for (int i = 0; i < nodesInRow; i++) {
    s1Nodes += [i];
  }

  for (int i = 0; i <= nodesInRow * (nodesInRow - 1); i += nodesInRow) {
    s1Nodes += [i];
  }

  for (int i = nodesInRow * (nodesInRow - 1);
      i < nodesInRow * nodesInRow;
      i++) {
    s1Nodes += [i];
  }

  for (int i = nodesInRow - 1; i < nodesInRow * nodesInRow; i += nodesInRow) {
    s1Nodes += [i];
  }

  return s1Nodes;
}

void main(List<String> args) async {
  var beg = 0.0;
  var end = 4.0;
  var step = 1.0;

  try {
    if (args.length == 1) {
      step = double.parse(args[0]);
    } else if (args.length == 3) {
      beg = double.parse(args[0]);
      end = double.parse(args[1]);
      step = double.parse(args[2]);
    }
  } catch (e) {
    print("Неправильное использование скрипта.");
    print("Правильное использование: ");
    print(
        "  - без параметров (в таком случае beg = 0.0, end = 4.0, step = 1.0)");
    print(
        "  - с одним параметром (в таком случае в качестве параметра передаётся step)");
    print(
        "  - с тремя параметрами (в таком случае вручную задаются все 3 параметра по порядку)");
    print("");
    print("Пример использования:");
    print(
        "  - ./generateGrid.exe 0.5             - задаст сетку от 0 до 4 с шагом 0.5");
    print(
        "  - ./generateGrid.exe 0.0 2.0 0.5     - задаст сетку от 0 до 2 с шагом 0.5");
  }
  var (nodesInRow, nodes) = CreateNodes(beg, end, step);
  var finites = GenerateFinites(nodesInRow);
  var treugs = GenerateTreugs(finites);
  var s1Nodes = GenerateS1Nodes(nodesInRow);
  print("Число узлов: ${nodes.length}, число КЭ: ${treugs.length}");

  // Print nodes to file ./nodes.txt
  var out = "";
  out += nodes.length.toString() + "\n";
  for (var elem in nodes) {
    out += elem.toFileString() + "\n";
  }
  File("./point.txt")
    ..createSync()
    ..writeAsStringSync(out);

  out = "";
  out += treugs.length.toString() + "\n";
  for (var elem in treugs) {
    out += elem.toFileString() + "\n";
  }
  File("./triangles.txt")
    ..createSync()
    ..writeAsStringSync(out);

  out = "";
  out += s1Nodes.length.toString() + "\n";
  for (var elem in s1Nodes) {
    out += elem.toString().padLeft(3) + "   0\n";
  }
  File("./kraev1.txt")
    ..createSync()
    ..writeAsStringSync(out);
}

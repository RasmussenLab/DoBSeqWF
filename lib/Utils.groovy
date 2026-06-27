import nextflow.Channel
import static nextflow.Nextflow.*

class Utils {
    /**
     * Return a value channel carrying (mainFile, indexOrEmpty).
     * indexExt is a single extension like '.bai' or '.tbi'.
     */
    static def indexedFileChannel(String filePath, String indexExt = '.tbi') {
        if (!filePath) return Channel.value(tuple([], []))

        def main  = file(filePath)
        def index = file("${filePath}${indexExt}")

        return Channel.value(index.exists() ? tuple(main, index) : tuple(main, []))
    }
}